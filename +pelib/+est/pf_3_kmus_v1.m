%> PF_3_KMUS Particle Filter for performing sensor fusion on the trajectory of
%> three KMUs presumably worn on the body in the following configuration: mid
%> pelvis, left ankle, right ankle
%> In this state space model, the position and velocity of each kinematic
%> measurement unit (KMU) is estimated in 3D space by combining the
%> information from each KMU in a kalman filter. NOTE: pay special attention 
%> to units:
%> position (meters)
%> velocity (m/s)
%> acceleration (m/2^2)
%> uwb_mea (meters)
%>
%> Author: Luke Wicent Sy
%>
%> Inputs::
%>   fs - sampling frequency of the magnetic and inertial measurement units
%>   sigma_acc - user specified process noise, i.e., the standard deviation
%>               in the accelerometer measurements when subjected to a known
%>               acceleration
%>   x0        - the initial state in the GFR
%>   gfrAccMP - the acceleration of the mid-pelvis in the GFR
%>   gfrAccLA - the acceleration of the left ankle in the GFR
%>   gfrAccRA - the acceleration of the right ankle in the GFR
%>   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   qMP       - mid  pelvis orientation in the GFR (quaternion)
%>   qLA       - left  ankle orientation in the GFR (quaternion)
%>   qRA       - right ankle orientation in the GFR (quaternion)
%>   dPelvis   - pelvis width
%>   dRFemur   - right femur length
%>   dLFemur   - left femur length
%>   dRTibia   - right tibia length
%>   dLTibia   - left tibia length
%>   uwb_mea    - a structure containing the range measurements (m) between
%>   options   - struct containing the ff. settings:
function [ xhat_pri, xhat_pos, debug_dat ] = pf_3_kmus_v3(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)

    %% default configurations
    fOpt = struct('fs', 60, 'applyMeas', false, 'applyUwb', false, ...
        'applyAccBias', false, 'applyCstr', 0, ...
        'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 1e3, 'sigmaQOriLA', 1e3, 'sigmaQOriRA', 1e3, ...
        'sigmaROriMP', 1e-3, 'sigmaROriLA', 1e-3, 'sigmaROriRA', 1e-3, ...
        'sigmaRPosLA', 1e-2, 'sigmaRPosRA', 1e-2, 'sigmaRPosZMP', 1e-1, ...
        'sigmaRPosMPLimit', 1e2, 'sigmaRPosLALimit', 1e2, ...
        'sigmaRPosRALimit', 1e2, 'sigmaRPosLARecenter', 1e-3, ...
        'sigmaRVelMPLARA', 1e1, 'sigmaRPosMPLARA', 1e1, ...
        'sigmaCPos', 1e-2, ...
        'sigmaUwbMPLA', 1e1, 'sigmaUwbMPRA', 1e1, 'sigmaUwbLARA', 1e1, ...
        'sigmaUwbLLeg', 1e0, 'sigmaUwbRLeg', 1e0, ...
        'sigmaZuptMP', 1e-1, 'sigmaZuptLA', 1e-1, 'sigmaZuptRA', 1e-1, ...
        'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
        'alphaRKmin', 0, 'alphaRKmax', pi*8/9, ...
        'optimOptimalityTolerance', 1e-2, ...
        'optimConstraintTolerance', 1e-2, ...
        'optimMaxFunctionEvaluations', 1500, 'optimUseParallel', false, ...
        'sckfAlpha', 0.1, 'sckfThreshold', 100, 'sckfMaxIter', 500);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% initialization
    idxPosMP = 1:3;
    idxVelMP = 4:6;
    idxPosVelMP = [idxPosMP idxVelMP];
	idxOriMP = 7:10;
    idxOriLT = 11:13;
    idxOriRT = 14:16;
    idxOriLK = 17:17;
    idxOriRK = 18:18;  
    nStates = 18;
    debug_dat = {};
    
    % initialise state vector (must be column)
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', nStates});
    fs = fOpt.fs;
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt*dt;        % local variable for readability
    I_N = eye(nStates);
    
    % state transition matrix encodes the relationship between previous state
    % estimate and current state estimate
    G = zeros(nStates, 3);
    G(idxPosMP, 1:3) = dt2.*eye(3);
    G(idxVelMP, 1:3) = dt .*eye(3);

    sigmaQAccMP = repelem((fOpt.sigmaQAccMP)^2, 3);
    % Initialise process noise covariance
    if islogical(P0) && ~P0
        Q = G * diag(sigmaQAccMP) * G';   
        P0 = Q;
    elseif isscalar(P0)
        P0 = P0*I_N;
    end
    
    nMeasure = 12;
   
    % check that all accelerometer measurements are equal dimensions
    [nSamples, ~] = size(gfrAccMP);
    validateattributes(gfrAccMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(qMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    
    % local variable assignment for readability
    u_k = [gfrAccMP, gfrAccLA, gfrAccRA]';
    y_k = [qMP, qLA, qRA]';
    uwbMPLARA = [uwb_mea.left_tibia_mid_pelvis,...
                 uwb_mea.mid_pelvis_right_tibia,...
                 uwb_mea.left_tibia_right_tibia];
                 
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);
    
    pf = particleFilter(@stateTransitionFcn001, @measLikelihoodFcn001);
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    nParticles = 1000;
    initialize(pf, nParticles, x0, P0);

    for n = 1:nSamples
        % Predict next position. Resample particles if necessary.
        [xhat_pri(n,:), P_pri(:,:,n)] = predict(pf, fs, ...
                                gfrAccMP(n, :)', qMP(n, :)', Q);
        % Generate dot measurement with random noise. This is
        % equivalent to the observation step.
%         measurement(i,:) = dot(i,:) + noise*(rand([1 2])-noise/2);
        measurement = 0;
        % Correct position based on the given measurement to get best estimation.
        % Actual dot position is not used. Store corrected position in data array.
        [xhat_pos(n,:), P_pos(:,:,n)] = correct(pf, measurement);
    end
end

function predParticles = stateTransitionFcn001(prevParticles, ...
                                            fs, gfrAccMP, qMP, procNoise)
    [nStates, nParticles] = size(prevParticles);  
    idxPosMP = 1:3; idxVelMP = 4:6; idxOriMP = 7:10;
    idxPosVelMP = 1:6;

    dt = 1.0/fs; % [s] Sample time
    dt2 = 0.5*dt*dt;      % local variable for readability
    
    F = eye(6, 6);
    F(idxPosMP, idxVelMP) = dt.*eye(3);   
    G = zeros(6, 3);
    G(idxPosMP, :) = dt2.*eye(3);
    G(idxVelMP, :) = dt .*eye(3);

    predParticles = prevParticles;
    predParticles(idxPosVelMP, :) = F*prevParticles(1:6, :) + G*gfrAccMP;
    predParticles(idxOriMP, :) = repelem(qMP, 1, nParticles);

    predParticles = predParticles + procNoise * randn(nStates, nParticles);
end

function likelihood = measLikelihoodFcn001(predParticles, meas, varargin)
    numberOfMeasurements = 1; % Expected number of measurements

    % Validate the measurement
    validateattributes(meas, {'double'}, {'vector', 'numel', numberOfMeasurements}, ...
        'vdpMeasurementLikelihoodFcn', 'measurement');

    % Assume that measurements are subject to Gaussian distributed noise with
    % variance 0.016
    % Specify noise as covariance matrix
    measurementNoise = 0.016 * eye(numberOfMeasurements);

    % The measurement contains the first state variable. Get the first state of
    % all particles
    predictedMeasurement = predParticles(1,:);

    % Calculate error between predicted and actual measurement
    measurementError = bsxfun(@minus, predictedMeasurement, meas(:)');

    % Use measurement noise and take inner product
    measurementErrorProd = dot(measurementError, measurementNoise \ measurementError, 1);

    measurementErrorProd = zeros(1000, 1);
    % Convert error norms into likelihood measure. 
    % Evaluate the PDF of the multivariate normal distribution. A measurement
    % error of 0 results in the highest possible likelihood.
    likelihood = 1/sqrt((2*pi).^numberOfMeasurements * det(measurementNoise)) * exp(-0.5 * measurementErrorProd);
end