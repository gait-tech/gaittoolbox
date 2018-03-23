function [ xhat_pri, xhat_pos, debug_dat ] = kf_3_kmus_v3(x0, P0, ...
    gfr_acc_MP, bIsStatMP, q_MP, ...
    gfr_acc_LA, bIsStatLA, q_LA, ...
    gfr_acc_RA, bIsStatRA, q_RA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)
% KF_3_KMUS Kalman Filter for performing sensor fusion on the trajectory of
% three KMUs presumably worn on the body in the following configuration: mid
% pelvis, left ankle, right ankle
% In this state space model, the position and velocity of each kinematic
% measurement unit (KMU) is estimated in 3D space by combining the
% information from each KMU in a kalman filter. NOTE: pay special attention 
% to units:
% position (meters)
% velocity (m/s)
% acceleration (m/2^2)
% uwb_mea (meters)
%
% Author: Michael Del Rosario, Luke Wicent Sy
%
% Inputs::
%   fs - sampling frequency of the magnetic and inertial measurement units
%   sigma_acc - user specified process noise, i.e., the standard deviation
%               in the accelerometer measurements when subjected to a known
%               acceleration
%   x0        - the initial state in the GFR
%   gfr_acc_MP - the acceleration of the mid-pelvis in the GFR
%   gfr_acc_LA - the acceleration of the left ankle in the GFR
%   gfr_acc_RA - the acceleration of the right ankle in the GFR
%   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   q_MP       - mid  pelvis orientation in the GFR (quaternion)
%   q_LA       - left  ankle orientation in the GFR (quaternion)
%   q_RA       - right ankle orientation in the GFR (quaternion)
%   dPelvis   - pelvis width
%   dRFemur   - right femur length
%   dLFemur   - left femur length
%   dRTibia   - right tibia length
%   dLTibia   - left tibia length
%   uwb_mea    - a structure containing the range measurements (m) between
%   options   - struct containing the ff. settings:
%       applyZupt - turn on/off zero velocity update. boolean
%       applyUwb - turn on/off uwb measurement update. boolean
%       applyAccBias - turn on/off acc bias in the model. boolean
%       applyConst - turn on/off constraints.
%   

    fOpt = struct('fs', 60, 'applyZupt', false, 'applyUwb', false, ...
        'applyAccBias', false, 'applyConst', 0, ...
        'sigmaAccMP', 0.5, 'sigmaAccLA', 0.5, 'sigmaAccRA', 0.5, ...
        'sigmaOriMP', 1e-2, 'sigmaOriLA', 1e-2, 'sigmaOriRA', 1e-2, ...
        'sigmaUwbMPLA', 0.2, 'sigmaUwbMPRA', 0.2, 'sigmaUwbLARA', 0.1, ...
        'sigmaZuptMP', 0.5, 'sigmaZuptLA', 0.5, 'sigmaZuptRA', 0.5, ...
        'optimOptimalityTolerance', 1e-2, ...
        'optimConstraintTolerance', 1e-2, ...
        'optimMaxFunctionEvaluations', 1500);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10 % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    idxVelLA = 14:16; % column idx corresponding to the left ankle velocity
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    idxVelRA = 24:26; % column idx corresponding to the right ankle velocity
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation
    idxMOriMP = 1:4;
    idxMOriLA = 5:8;
    idxMOriRA = 9:12;
    
    nStates = 30;
  
    if fOpt.applyAccBias
        idxAccBiasMP = 31:33; % column idx corresponding to the mid-pelvis acc bias
        idxAccBiasLA = 34:36; % column idx corresponding to the left ankle acc bias
        idxAccBiasRA = 37:39; % column idx corresponding to the right ankle acc bias
        nStates = 39;
    end
    
    % initialise state vector (must be column)
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', nStates});
    xhat = x0;

    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt*dt;        % local variable for readability
    I_N = eye(nStates);
    
    % state transition matrix encodes the relationship between previous state
    % estimate and current state estimate
    F = eye(nStates,nStates);
    % x = x(t-1) + v(t-1)*dt + 0.5*a(t)*dt^2
    % x = A*x + B*u
    F(idxPosMP,idxVelMP) = dt.*eye(3); % mid pelvis
    F(idxPosLA,idxVelLA) = dt.*eye(3); % left ankle
    F(idxPosRA,idxVelRA) = dt.*eye(3); % right ankle

    if applyAccBias
        F(idxPosMP, idxAccBiasMP) = -dt2.*eye(3);
        F(idxPosLA, idxAccBiasLA) = -dt2.*eye(3);
        F(idxPosRA, idxAccBiasRA) = -dt2.*eye(3);
        F(idxVelMP, idxAccBiasMP) = -dt.*eye(3);
        F(idxVelLA, idxAccBiasLA) = -dt.*eye(3);
        F(idxVelRA, idxAccBiasRA) = -dt.*eye(3);
    end
    
    G = zeros(nStates,9);
    G(idxPosMP, 1:3) = dt2.*eye(3);
    G(idxVelMP, 1:3) = dt .*eye(3);
    G(idxPosLA, 4:6) = dt2.*eye(3);
    G(idxVelLA, 4:6) = dt .*eye(3);
    G(idxPosRA, 7:9) = dt2.*eye(3);
    G(idxVelRA, 7:9) = dt .*eye(3);

    % Initialise process noise covariance
    Q = diag(repelem([(fOpt.sigmaAccMP)^2 (fOpt.sigmaAccLA)^2 (fOpt.sigmaAccRA)^2], 3));
    Q = G * Q * G';

    % initialise covariance in the state estimate
    if islogical(P0) && ~P0
        P_plus = Q;
    else
        P_plus = P0;
    end
    
    nMeasure = 12;
    H = zeros(nMeasure, nStates);
    H(idxMOriMP, idxOriMP) = eye(4, 4);
    H(idxMOriLA, idxOriLA) = eye(4, 4);
    H(idxMOriRA, idxOriRA) = eye(4, 4);

    Rdiag = repelem([(fOpt.sigmaOriMP)^2 (fOpt.sigmaOriLA)^2 (fOpt.sigmaOriRA)^2], 4);
    R = diag(Rdiag);
    
    % check that all accelerometer measurements are equal dimensions
    [nSamples, ~] = size(gfr_acc_MP);
    validateattributes(gfr_acc_MP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfr_acc_LA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfr_acc_RA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(q_MP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(q_LA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(q_RA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    
    % local variable assignment for readability
    u_k = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA]';
    y_k = [q_MP, q_LA, q_RA]';
    
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);

    debug_dat = struct;
    debug_dat.LFEO = nan(nSamples, 3); debug_dat.RFEO = nan(nSamples, 3);
    debug_dat.LFEP = nan(nSamples, 3); debug_dat.RFEP = nan(nSamples, 3);
    debug_dat.qLFemur = nan(nSamples, 4); debug_dat.qRFemur = nan(nSamples, 4);
    
    debug_dat.predState = nan(nSamples, nStates);
    debug_dat.zuptStateL = nan(nSamples, nStates);
    debug_dat.zuptStateR = nan(nSamples, nStates);
    debug_dat.uwbuptState = nan(nSamples, nStates);
    debug_dat.cstrState = nan(nSamples, nStates);
    
    if fOpt.applyZupt
        idxMVelMP = nMeasure+1:nMeasure+3;
        idxMVelLA = nMeasure+4:nMeasure+6;
        idxMVelRA = nMeasure+7:nMeasure+9;
        nMeasure = nMeasure+9;
        
        H(end+1:end+9, nStates) = zeros(9, nStates);
        H(idxMVelMP, idxVelMP) = eye(3);
        H(idxMVelLA, idxVelLA) = eye(3);
        H(idxMVelRA, idxVelRA) = eye(3);
        
        Rdiag = diag(R);
        Rdiag(end+1:end+9) = repelem([(fOpt.sigmaZuptMP)^2 ...
            (fOpt.sigmaZuptLA)^2 (fOpt.sigmaZuptRA)^2], 3);
        R = diag(Rdiag);
        
        y_k(end+1:end+9, 1:nStates) = zeros(9, nSamples);
    end
    
    if fOpt.applyUwb
        % pseudo UWB measurements corresponding to the euclidean distance between
        % pairs of KMUs.
        % NOTE: the column order of these measurements is important. Assume the
        % following order unless stated otherwise:
        %  uwb_MP_LA_RA =  ['mid-pelvis to left ankle',
        %                   'mid-pelvis to right ankle',
        %                   'left ankle to right ankle'];
        uwb_MP_LA_RA = [uwb_mea.left_tibia_midPelvis,...
                        uwb_mea.midPelvis_right_tibia,...
                        uwb_mea.left_tibia_right_tibia];

        % specify the measurement noise in the UWB measurements, these may be
        % different for each KMU sensor pair. It is likely that the range between
        % feet/ankles will be the most accurate due to less "no line of sight"
        % periods. Note: units on sigma_uwb = meters
        R_uwb = zeros(3,3);
        R_uwb(1,1) = (fOpt.sigmaUwbMPLA)^2;
        R_uwb(2,2) = (fOpt.sigmaUwbMPRA)^2;
        R_uwb(3,3) = (fOpt.sigmaUwbLARA)^2;
    end
    
    if fOpt.applyConst
        D = [-eye(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
             -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3)];
     
        optimOpt = optimoptions('fmincon', 'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'OptimalityTolerance', fOpt.optimOptimalityTolerance, ...
            'ConstraintTolerance', fOpt.optimConstraintTolerance, ...
            'MaxFunctionEvaluations', fOpt.optimMaxFunctionEvaluations);
    end

    for n = 1:nSamples

    %% -----------------------------------------------------------------------
    % Prediction Step using accelerometer measurements ----    
%         if fOpt.applyAccBias
%             F(idxPosMP,idxAccBiasMP) = -dt2.*PELV_CS;
%             F(idxPosLA,idxAccBiasLA) = -dt2.*LTIB_CS;
%             F(idxPosRA,idxAccBiasRA) = -dt2.*RTIB_CS;
%             F(idxVelMP,idxAccBiasMP) = -dt.*PELV_CS;
%             F(idxVelLA,idxAccBiasLA) = -dt.*LTIB_CS;
%             F(idxVelRA,idxAccBiasRA) = -dt.*RTIB_CS;
%         end
        
        xhat = F * xhat + G * u_k(:,n) ;
        P_min= F * P_plus * F' + Q;
        xhat_pri(n,:) = xhat;
        P_pri(:,:,n)  = P_min;
        
        debug_dat.predState(n,:) = xhat;

    %% ------------------------------------------------------------------------
    % Measurement update step
    % matrices beginnning with 'H_' are the 'observation matrices' that map
    % the variables in the state estimate vector, xhat, to the measurement
    % domain. In this case we are using
        idx = [idxMOriMP idxMOriLA idxMOriRA];
        if fOpt.applyZupt
            ctrZUPT = 0;
            if bIsStatMP(n)
                ctrZUPT = ctrZUPT+1;
                idx(end+1:end+3) = idxMVelMP;
            end
            if bIsStatLA(n)
                ctrZUPT = ctrZUPT+1;
                idx(end+1:end+3) = idxMVelLA;
            end
            if bIsStatRA(n)
                ctrZUPT = ctrZUPT+1;
                idx(end+1:end+3) = idxMVelRA;
            end
        end
        res = y_k(idx, n) - H(idx, :) * xhat;
        K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx))

        xhat = xhat + K * y_k(idx, n);
        P_min1 = (I_N - K * H(idx, :)) * P_min;

        if fOpt.applyZupt
            if bIsStatLA(n)
                debug_dat.zuptStateL(n,:) = xhat;
            end
            if bIsStatRA(n)
                debug_dat.zuptStateR(n,:) = xhat;
            end
        end
        
    %% ---- Kalman Filter Update Step using UWB measurements ---- 
    % this correction step should be done last
        if fOpt.applyUwb
            diff_MP_LA = xhat(idxPosMP)' - xhat(idxPosLA)';
            diff_MP_RA = xhat(idxPosMP)' - xhat(idxPosRA)';
            diff_LA_RA = xhat(idxPosLA)' - xhat(idxPosRA)';
            % the observation model
            h_uwb_est = [vecnormalize( diff_MP_LA );
                         vecnormalize( diff_MP_RA );
                         vecnormalize( diff_LA_RA );];
            % calculate the measurement residual, i.e., the difference between the
            % predicted measurements from the observation model, h(xhat), and the
            % measurements obtained directly from the UWB sensors.
            y_innovation = uwb_MP_LA_RA(n,:)' - h_uwb_est;

            H_uwb = zeros(3, nStates);
            % Jacobian of observation model with respect to elements of the state
            % estimate xhat, the order of these matrix elements MUST be preserved
            H_uwb(1,1) = diff_MP_LA(1)/h_uwb_est(1);
            H_uwb(1,2) = diff_MP_LA(2)/h_uwb_est(1);
            H_uwb(1,3) = diff_MP_LA(3)/h_uwb_est(1);
            H_uwb(1,7) = -diff_MP_LA(1)/h_uwb_est(1);
            H_uwb(1,8) = -diff_MP_LA(2)/h_uwb_est(1);
            H_uwb(1,9) = -diff_MP_LA(3)/h_uwb_est(1);    

            H_uwb(2,1)  =  diff_MP_RA(1)/h_uwb_est(2);
            H_uwb(2,2)  =  diff_MP_RA(2)/h_uwb_est(2);
            H_uwb(2,3)  =  diff_MP_RA(3)/h_uwb_est(2);
            H_uwb(2,13) = -diff_MP_RA(1)/h_uwb_est(2);
            H_uwb(2,14) = -diff_MP_RA(2)/h_uwb_est(2);
            H_uwb(2,15) = -diff_MP_RA(3)/h_uwb_est(2);

            H_uwb(3,7)  =  diff_LA_RA(1)/h_uwb_est(3);
            H_uwb(3,8)  =  diff_LA_RA(2)/h_uwb_est(3);
            H_uwb(3,9)  =  diff_LA_RA(3)/h_uwb_est(3);
            H_uwb(3,13) = -diff_LA_RA(1)/h_uwb_est(3);
            H_uwb(3,14) = -diff_LA_RA(2)/h_uwb_est(3);
            H_uwb(3,15) = -diff_LA_RA(3)/h_uwb_est(3);    
            % Calculate the covariance in the measurement residual
            S    = ((H_uwb * P_min1) * H_uwb') + R_uwb; % scalar
            % Calculate the kalman gain
            K    = P_min1 * H_uwb' * S^(-1);
            % Update state estimate
            xhat = xhat + K * y_innovation;
            % Update Covariance in the state estimate
            P_plus    = (I_N - K * H) * P_min1;

            debug_dat.uwbuptState(n,:) = xhat;
        else
            P_plus = P_min1;
        end

    %% -----------------------------------------------------------------------
    % Constraint update step ---- 
        LTIB_CS = quat2rotm(xhat(idxOriLA,1));
        RTIB_CS = quat2rotm(xhat(idxOriRA,1));
        PELV_CS = quat2rotm(xhat(idxOriMP,1));
        % Test frankenstein constraint
        if fOpt.applyConst == 1
            % calculate the location of the knee
            LKNE = xhat(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
            RKNE = xhat(idxPosRA,1) + dRTibia*RTIB_CS(:,3);

            % calculate the z axis of the femur
            LFEM_z = xhat(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
            RFEM_z = xhat(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;

            % calculate the z axis of the tibia
            LTIB_z = LTIB_CS(:,3);
            RTIB_z = RTIB_CS(:,3);

            % calculate alpha_lk and alpha_rk
            alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
            alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));

            % setup the constraint equations
            d_k = [ (dPelvis/2*PELV_CS(:,2) ...
                     -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
                     +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
                     -dLTibia*LTIB_CS(:,3)) ; ...
                    (-dPelvis/2*PELV_CS(:,2)+ ...
                     -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
                     +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
                     -dRTibia*RTIB_CS(:,3)) ];
            
            Kk = D'*(D*D')^(-1);
            res = d_k - D * xhat;
            dx = Kk*(res);
            xhat = xhat + dx;
            
            debug_dat.cstrStateRes(n,:) = res;
            debug_dat.cstrState(n,:) = xhat;
            debug_dat.cstrStateKk(:,:,n) = Kk;
        end

        debug_dat.LFEO(n, :) = xhat(idxPosLA) + dLTibia * LTIB_CS(:, 3);
        debug_dat.RFEO(n, :) = xhat(idxPosRA) + dRTibia * RTIB_CS(:, 3);
        debug_dat.LFEP(n, :) = xhat(idxPosMP) + dPelvis/2 * PELV_CS(:, 2);
        debug_dat.RFEP(n, :) = xhat(idxPosMP) - dPelvis/2 * PELV_CS(:, 2);
        LFEM_z = (debug_dat.LFEP(n,:)-debug_dat.LFEO(n,:))';
        LFEM_y = LTIB_CS(:,2);
        LFEM_x = cross(LFEM_y, LFEM_z);
        RFEM_z = (debug_dat.RFEP(n,:)-debug_dat.RFEO(n,:))';
        RFEM_y = RTIB_CS(:,2);
        RFEM_x = cross(RFEM_y, RFEM_z);
        debug_dat.qLTH(n, :) = rotm2quat([LFEM_x LFEM_y LFEM_z]);
        debug_dat.qRTH(n, :) = rotm2quat([RFEM_x RFEM_y RFEM_z]);

        xhat_pos(n, :) = xhat;
        P_pos(:, :, n)  = P_plus;
    end
end

function y = L2Dist(x, x0, S)
%x^2 is monotomically increasing at any point not at 0, so traditional
%L2 norm involving sqrt is unnecessary, can use x^2 to find same
%location of min cost in constrained region with less computational
%cost.
%using inverse of covariance rather than I will make state est over time
%less smooth but more accurate over the average of the interval

    res = (x-x0);
    res = S\res; %add res*inv(S) to scale cost by certainty
    n = 1000; %have also tried 2,4,6,8,14,100,1000 around >= 14 greatly increases speed of finding solution, not much difference etween 100 and 1000
    y = (res'*res)^n;
end

function [c, ceq] = hjc_kneeaug_nonlcon(x, Rpelv2, dPelvis, ...
    dLFemur, dRFemur, dLTibia, dRTibia)
    LFEM_z = x(1:3,1)+dPelvis/2*Rpelv2-x(7:9,1);
    RFEM_z = x(1:3,1)-dPelvis/2*Rpelv2-x(13:15,1);
    LTIB_z = x(7:9,1)-x(19:21,1);
    RTIB_z = x(13:15,1)-x(25:27,1);

    ceq = [LFEM_z'*LFEM_z - dLFemur^2;
           RFEM_z'*RFEM_z - dRFemur^2;
           LTIB_z'*LTIB_z - dLTibia^2;
           RTIB_z'*RTIB_z - dRTibia^2];
    c = [];
end