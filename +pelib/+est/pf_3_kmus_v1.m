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
function [ xhat_pri, xhat_pos, debug_dat ] = pf_3_kmus_v1(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)

    %% default configurations
    fOpt = struct('fs', 60, 'applyMeas', 1, 'applyUwb', false, ...
        'applyAccBias', false, 'applyPred', 1, ...
        'sigmaQAccMP', 1, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 1e3, 'sigmaQOriLA', 1e3, 'sigmaQOriRA', 1e3, ...
        'sigmaROriMP', 1e-3, 'sigmaROriLA', 1e-3, 'sigmaROriRA', 1e-3, ...
        'sigmaROriLSK', 1e-2, 'sigmaROriRSK', 1e-2, ...
        'sigmaRAccLA', 1e1, 'sigmaRAccRA', 1e1, ...
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
        P0 = G * diag(sigmaQAccMP) * G';
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
    % gfrAccMP, gfrAccLA, gfrAccRA;
    y_k = [qMP, qLA, qRA]';
    uwbMPLARA = [uwb_mea.left_tibia_mid_pelvis,...
                 uwb_mea.mid_pelvis_right_tibia,...
                 uwb_mea.left_tibia_right_tibia];
    dBody = struct('RPV', dPelvis, 'LTH', dLFemur, 'RTH', dRFemur, ...
                   'LSK', dLTibia, 'RSK', dRTibia);
    body0 = pelib.grBody.generateBodyFromJointAngles(x0(1:3,1)', ...
            x0(7:10,1)', x0(11:13,1)', x0(14:16,1)', ...
            [0 x0(17,1) 0], [0 x0(18,1) 0], ...
            dBody.RPV, dBody.LTH, dBody.RTH, dBody.LSK, dBody.RSK);
    zFloor = mean([body0.LTIO(1,3), body0.RTIO(1,3)]);
    
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);
    
    applyPred = fOpt.applyPred;
    applyMeas = fOpt.applyMeas;
    switch fOpt.applyPred
        case 1
            stateTransitionFcn = @stateTransitionFcn001;
            Q = G * diag(sigmaQAccMP) * G' + diag([zeros(10, 1); deg2rad(1)*ones(8,1)]);
        case 2
            stateTransitionFcn = @stateTransitionFcn002;
            sigma = struct('accMP', fOpt.sigmaQAccMP, ...
                           'oriLK', deg2rad(0.1), 'oriRK', deg2rad(0.1));
    end
    
    if (applyMeas >= 1 && applyMeas <= 7)
        measLikelihoodFcn = @measLikelihoodFcn001;
%         R = diag([1e-2*ones(14,1); (fOpt.sigmaQAccLA)^2*ones(3,1); ...
%               (fOpt.sigmaQAccRA)^2*ones(3,1)]);
        R = diag([(fOpt.sigmaROriLSK).^2*ones(3,1); ...
                  (fOpt.sigmaROriRSK).^2*ones(3,1);
                  (fOpt.sigmaZuptLA).^2*ones(3,1); (fOpt.sigmaRPosLA).^2; ...
                  (fOpt.sigmaZuptRA).^2*ones(3,1); (fOpt.sigmaRPosRA).^2; ...
                  (fOpt.sigmaRAccLA).^2*ones(3,1)
                  (fOpt.sigmaRAccRA).^2*ones(3,1)]);
    end
    
    pf = particleFilter(stateTransitionFcn, measLikelihoodFcn);
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    nParticles = 1000;
    initialize(pf, nParticles, x0, P0);

    for n = 1:nSamples
        % Predict next position. Resample particles if necessary.
        switch fOpt.applyPred
            case 1
                [xhat_pri(n,:), P_pri(:,:,n)] = predict(pf, fs, Q, ...
                                        gfrAccMP(n, :)', qMP(n, :)');
            case 2
                qOri = struct('qRPV', qMP(n, :), ...
                              'qLSK', qLA(n, :), 'qRSK', qRA(n, :));
                [xhat_pri(n,:), P_pri(:,:,n)] = predict(pf, fs, sigma, ...
                                        gfrAccMP(n, :)', qOri);
        end
        % Generate dot measurement with random noise. This is
        % equivalent to the observation step.
%         measurement(i,:) = dot(i,:) + noise*(rand([1 2])-noise/2);
        % Correct position based on the given measurement to get best estimation.
        % Actual dot position is not used. Store corrected position in data array.
        if (applyMeas >= 1 && applyMeas <= 7)
%             gfrAcc = struct('MP', gfrAccMP(n, :), 'LA', gfrAccLA(n, :), ...
%                         'RA', gfrAccRA(n, :));
            step = [bIsStatLA(n, 1) bIsStatRA(n, 1)];
            if n-1 >= 1, past1Particles = xhat_pos(n-1,:)';
            else, past1Particles = x0; end
            if n-2 >= 1, past2Particles = xhat_pos(n-2,:)';
            else, past2Particles = x0; end
        
            meas = [qLA(n, :) qRA(n, :) gfrAccLA(n, :) gfrAccRA(n,:)];
            [xhat_pos(n,:), P_pos(:,:,n)] = correct(pf, meas, fs, R, applyMeas, ...
                    dBody, step, zFloor, past1Particles, past2Particles);
        end
    end
end

function predParticles = stateTransitionFcn001(prevParticles, ...
                                            fs, Q, gfrAccMP, qMP)
    % starting point: MP pos, vel and ori
    % varies: MP pos, vel; Hip and knee joint angles
    
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

    predParticles = predParticles + Q * randn(nStates, nParticles);
    
    %% knee constraint (no hyperextension)
    predParticles(17:18, :) = max(predParticles(17:18, :), 0);
end

function predParticles = stateTransitionFcn002(prevParticles, ...
                                            fs, sigma, gfrAccMP, qOri)
    % starting point: MP pos, vel; RPV LSK RSK ori
    % varies: MP pos, vel; Knee joint angles
    
    [nStates, nParticles] = size(prevParticles);  
    idxPosMP = 1:3; idxVelMP = 4:6; idxOriRPV = 7:10;
    idxPosVelMP = 1:6; idxOriLT = 11:13; idxOriRT = 14:16;
    idxOriLK = 17; idxOriRK = 18; idxOriLRK = 17:18;
    seq = [2 1 3];
    
    dt = 1.0/fs; % [s] Sample time
    dt2 = 0.5*dt*dt;      % local variable for readability
    
    F = eye(6, 6);
    F(idxPosMP, idxVelMP) = dt.*eye(3);   
    G = zeros(6, 3);
    G(idxPosMP, :) = dt2.*eye(3);
    G(idxVelMP, :) = dt .*eye(3);

    predParticles = prevParticles;
    predParticles(idxPosVelMP, :) = F*prevParticles(1:6, :) + ...
                G*(gfrAccMP + sigma.accMP*randn(3, nParticles));
    predParticles(idxOriRPV, :) = repelem(qOri.qRPV', 1, nParticles);
    predParticles(idxOriLK, :) = predParticles(idxOriLK, :) + ...
                + sigma.oriLK*randn(1, nParticles);
    predParticles(idxOriRK, :) = predParticles(idxOriRK, :) + ...
                + sigma.oriRK*randn(1, nParticles);        
    % knee constraint (no hyperextension)
    predParticles(idxOriLRK, :) = max(predParticles(17:18, :), 0);
    
    % solve hip angles from knee angle and tibia orientation
    zN = zeros(nParticles, 1);
    LKneeAngle = [zN, predParticles(idxOriLK, :)', zN].*[-1 1 -1];
    RKneeAngle = [zN, predParticles(idxOriRK, :)', zN];
    qLTH  = pelib.grBody.calcProxRotm(qOri.qLSK, LKneeAngle(:,seq));
    qRTH  = pelib.grBody.calcProxRotm(qOri.qRSK, RKneeAngle(:,seq));
    LHipAngles = pelib.grBody.calcJointAngles(predParticles(idxOriRPV, :)', qLTH);
    RHipAngles = pelib.grBody.calcJointAngles(predParticles(idxOriRPV, :)', qRTH);
    predParticles(idxOriLT,:) = (LHipAngles(:,seq).*[-1 -1 -1])';
    predParticles(idxOriRT,:) = (RHipAngles(:,seq).*[ 1 -1  1])';    
end

function likelihood = measLikelihoodFcn001(predParticles, meas, fs, R, ...
        applyMeas, d, step, zFloor, past1Particles, past2Particles)
    %% initialization
    bodyList = {};
    [nStates, nParticles] = size(predParticles);  
    particleList = {predParticles, past1Particles, past2Particles};
    
    for i=1:3
        p = particleList{i}';
        zN = zeros(size(p, 1), 1);
        
        bodyList{i} = pelib.grBody.generateBodyFromJointAngles(p(:,1:3), ...
            p(:,7:10), p(:,11:13), p(:,14:16), [zN p(:,17) zN], [zN p(:,18) zN], ...
            d.RPV, d.LTH, d.RTH, d.LSK, d.RSK);
    end
    seq = 'YXZ';
    dt = 1.0/fs; dt2 = dt^2;
    idxMOri = 1:8; idxROri = 1:6;
    idxRLStep = 7:10;
    idxRRStep = 11:14;
    idxMAccLA = 9:11; idxRAccLA = 15:17;
    idxMAccRA = 12:14; idxRAccRA = 18:20;
    modTenApplyMeas = mod(applyMeas, 10);
    
    %% LTIB and RTIB orientation
    nMeas = 0;
    if bitand(modTenApplyMeas, 1)
        nMeas = nMeas + 6; % Expected number of measurements
    end
    if bitand(modTenApplyMeas, 2)
        if step(1), nMeas = nMeas + 4; end % left  step detection
        if step(2), nMeas = nMeas + 4; end % right step detection
    end
    if bitand(modTenApplyMeas, 4)
        nMeas = nMeas + 6; 
    end
    
    idxR = zeros(nMeas, 1);
    tempIdxR = 1;
    
    %% qLSK and qRSK
    if bitand(modTenApplyMeas, 1)
        idxR(tempIdxR:tempIdxR+5) = idxROri;
        tempIdxR = tempIdxR+6;
        [d1, d2, d3] = quat2angle(quatmultiply(bodyList{1}.qLSK, quatconj(meas(:,idxMOri(1:4)))), seq);
        [d4, d5, d6] = quat2angle(quatmultiply(bodyList{1}.qRSK, quatconj(meas(:,idxMOri(5:8)))), seq);
        diffOri = [d1, d2, d3, d4, d5, d6]';
    else
        diffOri = [];
    end
    
    %% ZUPT and floor assumption
    if bitand(modTenApplyMeas, 2) && step(1)
        diffLStep = [(bodyList{1}.LTIO-bodyList{2}.LTIO)'./dt; ...
                     bodyList{1}.LTIO(:,3)'-zFloor];
        idxR(tempIdxR:tempIdxR+3) = idxRLStep;
        tempIdxR = tempIdxR+4;
    else
        diffLStep = [];
    end
    if bitand(modTenApplyMeas, 2) && step(2)
        diffRStep = [(bodyList{1}.RTIO-bodyList{2}.RTIO)'./dt; ...
                     bodyList{1}.RTIO(:,3)'-zFloor];
        idxR(tempIdxR:tempIdxR+3) = idxRRStep;
        tempIdxR = tempIdxR+4;
    else
        diffRStep = []; 
    end

    %% Acceleration
    % adding acceleration breaks the already non converging solution
    if bitand(modTenApplyMeas, 4)
        predAccLA = (bodyList{1}.LTIO-2*bodyList{2}.LTIO+bodyList{3}.LTIO)/dt2;
        predAccRA = (bodyList{1}.RTIO-2*bodyList{2}.RTIO+bodyList{3}.RTIO)/dt2;
        diffAcc = [(predAccLA - meas(:,idxMAccLA)) (predAccRA - meas(:,idxMAccRA))]';
        idxR(tempIdxR:tempIdxR+2) = idxRAccLA;
        idxR(tempIdxR+3:tempIdxR+5) = idxRAccRA;
        tempIdxR = tempIdxR+6;
    else
        diffAcc = [];
    end
    
    %% culmination
    measError = [diffOri; diffLStep; diffRStep; diffAcc];
    R2 = R(idxR,idxR);
    % Use measurement noise and take inner product
    measErrorProd = dot(measError, R2 \ measError, 1);

    % Convert error norms into likelihood measure. 
    % Evaluate the PDF of the multivariate normal distribution. A measurement
    % error of 0 results in the highest possible likelihood.
    if size(measError, 2) == 0
        likelihood = zeros(1, nParticles);
    else
        likelihood  = 1/sqrt((2*pi).^nMeas * det(R2)) * exp(-0.5 * measErrorProd);
    end
end