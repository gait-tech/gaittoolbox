function [ xhat_pri, xhat_con, debug_dat ] = kf_3_kmus_v3(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
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
%   gfrAccMP - the acceleration of the mid-pelvis in the GFR
%   gfrAccLA - the acceleration of the left ankle in the GFR
%   gfrAccRA - the acceleration of the right ankle in the GFR
%   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   qMP       - mid  pelvis orientation in the GFR (quaternion)
%   qLA       - left  ankle orientation in the GFR (quaternion)
%   qRA       - right ankle orientation in the GFR (quaternion)
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
%           001: projection (W=I) assuming perfect orientation
%           101: fmincon

    fOpt = struct('fs', 60, 'applyZupt', false, 'applyUwb', false, ...
        'applyAccBias', false, 'applyConst', 0, ...
        'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 1e5, 'sigmaQOriLA', 1e5, 'sigmaQOriRA', 1e5, ...
        'sigmaROriMP', 1e-1, 'sigmaROriLA', 1e-1, 'sigmaROriRA', 1e-1, ...
        'sigmaUwbMPLA', 0.2, 'sigmaUwbMPRA', 0.2, 'sigmaUwbLARA', 0.1, ...
        'sigmaZuptMP', 1e-4, 'sigmaZuptLA', 1e-4, 'sigmaZuptRA', 1e-4, ...
        'optimOptimalityTolerance', 1e-2, ...
        'optimConstraintTolerance', 1e-2, ...
        'optimMaxFunctionEvaluations', 1500, 'optimUseParallel', false);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
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
    x_tilde = x0;

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

    if fOpt.applyAccBias
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
    Q = diag(repelem([(fOpt.sigmaQAccMP)^2 (fOpt.sigmaQAccLA)^2 (fOpt.sigmaQAccRA)^2], 3));
    Qori = diag([zeros(1,6) repelem((fOpt.sigmaQOriMP)^2, 1, 4) ...
                 zeros(1,6) repelem((fOpt.sigmaQOriLA)^2, 1, 4) ...
                 zeros(1,6) repelem((fOpt.sigmaQOriRA)^2, 1, 4)]);
    Q = G * Q * G' + Qori;
    % initialise covariance in the state estimate
    if islogical(P0) && ~P0
        P_tilde = Q;
    elseif isscalar(P0)
        P_tilde = P0*I_N;
    else
        P_tilde = P0;
    end
    
    nMeasure = 12;
    H = zeros(nMeasure, nStates);
    H(idxMOriMP, idxOriMP) = eye(4, 4);
    H(idxMOriLA, idxOriLA) = eye(4, 4);
    H(idxMOriRA, idxOriRA) = eye(4, 4);

    Rdiag = repelem([(fOpt.sigmaROriMP)^2 (fOpt.sigmaROriLA)^2 (fOpt.sigmaROriRA)^2], 4);
    R = diag(Rdiag);
    
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
    
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);
    
    xhat_con = nan(nSamples, nStates);
    P_con    = nan(nStates, nStates, nSamples);

    debug_dat = struct;
    debug_dat.LFEO = nan(nSamples, 3); debug_dat.RFEO = nan(nSamples, 3);
    debug_dat.LFEP = nan(nSamples, 3); debug_dat.RFEP = nan(nSamples, 3);
    debug_dat.qLFemur = nan(nSamples, 4); debug_dat.qRFemur = nan(nSamples, 4);
    
    debug_dat.predState = nan(nSamples, nStates);
    debug_dat.predP = nan(nStates, nStates, nSamples);
%     debug_dat.zuptStateL = false(nSamples, 1);
%     debug_dat.zuptStateR = false(nSamples, 1);
    debug_dat.zuptState = nan(nSamples, nStates);
    debug_dat.zuptP = nan(nStates, nStates, nSamples);
    debug_dat.uwbuptState = nan(nSamples, nStates);
    debug_dat.uwbuptP = nan(nStates, nStates, nSamples);
    debug_dat.cstrState = nan(nSamples, nStates);
    debug_dat.cstrP = nan(nStates, nStates, nSamples);
    
    debug_dat.zuptStateL = bIsStatLA;
    debug_dat.zuptStateR = bIsStatRA;
    
    if fOpt.applyZupt
        idxMVelMP = nMeasure+1:nMeasure+3;
        idxMVelLA = nMeasure+4:nMeasure+6;
        idxMVelRA = nMeasure+7:nMeasure+9;
        nMeasure = nMeasure+9;
        
        H(end+1:end+9, :) = zeros(9, nStates);
        H(idxMVelMP, idxVelMP) = eye(3);
        H(idxMVelLA, idxVelLA) = eye(3);
        H(idxMVelRA, idxVelRA) = eye(3);
        
        Rdiag = diag(R);
        Rdiag(end+1:end+9) = repelem([(fOpt.sigmaZuptMP)^2 ...
            (fOpt.sigmaZuptLA)^2 (fOpt.sigmaZuptRA)^2], 3);
        R = diag(Rdiag);
        
        y_k(end+1:end+9, :) = zeros(9, nSamples);
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
    
    if fOpt.applyConst >= 1 && fOpt.applyConst <= 2
        D = zeros(6, nStates);
        D(1:3,idxPosMP) = -eye(3, 3);
        D(1:3,idxPosLA) = eye(3, 3);
        D(4:6,idxPosMP) = -eye(3, 3);
        D(4:6,idxPosRA) = eye(3, 3);
    elseif fOpt.applyConst >= 101 && fOpt.applyConst <= 106
        optimOpt = optimoptions('fmincon', 'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'OptimalityTolerance', fOpt.optimOptimalityTolerance, ...
            'ConstraintTolerance', fOpt.optimConstraintTolerance, ...
            'MaxFunctionEvaluations', fOpt.optimMaxFunctionEvaluations, ...
            'UseParallel', fOpt.optimUseParallel);
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
        
        x_min = F * x_tilde + G * u_k(:,n) ;
        P_min = F * P_tilde * F' + Q;
        xhat_pri(n,:) = x_min;
        P_pri(:,:,n)  = P_min;
        
        debug_dat.predState(n,:) = x_min;
        debug_dat.predP(:,:,n) = P_min;
        
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
        
        res = y_k(idx, n) - H(idx, :) * x_min;
        K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx));

        x_min1 = x_min + K * res;
        P_min1 = (I_N - K * H(idx, :)) * P_min;
        
        if fOpt.applyZupt
            debug_dat.zuptState(n,:) = x_min1;
            debug_dat.zuptP(:,:,n) = P_min1;
        end
        
    %% ---- Kalman Filter Update Step using UWB measurements ---- 
    % this correction step should be done last
        if fOpt.applyUwb
            diff_MP_LA = x_min1(idxPosMP)' - x_min1(idxPosLA)';
            diff_MP_RA = x_min1(idxPosMP)' - x_min1(idxPosRA)';
            diff_LA_RA = x_min1(idxPosLA)' - x_min1(idxPosRA)';
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
            x_plus = x_min1 + K * y_innovation;
            % Update Covariance in the state estimate
            P_plus    = (I_N - K * H) * P_min1;

            debug_dat.uwbuptState(n,:) = x_plus;
            debug_dat.uwbuptP(:,:,n) = P_plus;
        else
            x_plus = x_min1;
            P_plus = P_min1;
        end
        
        xhat_pos(n, :) = x_plus;
        P_pos(:, :, n)  = P_plus;
        
    %% -----------------------------------------------------------------------
    % Constraint update step ---- 
        LTIB_CS = quat2rotm(x_plus(idxOriLA,1)');
        RTIB_CS = quat2rotm(x_plus(idxOriRA,1)');
        PELV_CS = quat2rotm(x_plus(idxOriMP,1)');
        % Test frankenstein constraint
        
        if fOpt.applyConst == 0
            x_tilde = x_plus;
            P_tilde = P_plus;
        elseif fOpt.applyConst >= 1 && fOpt.applyConst <= 2 % projection (W=I) assuming perfect orientation
            % calculate the location of the knee
            LKNE = x_plus(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
            RKNE = x_plus(idxPosRA,1) + dRTibia*RTIB_CS(:,3);

            % calculate the z axis of the femur
            LFEM_z = x_plus(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
            RFEM_z = x_plus(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;

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
            
            switch (fOpt.applyConst)
                case 1 % least squares estimate + no P update
                    Kk = D'*(D*D')^(-1);
                case 2 % maximum probability estimate + no P update
                    Kk = P_plus*D'*(D*P_plus*D')^(-1);
            end
            
            res = d_k - D * x_plus;
            dx = Kk*(res);
            x_tilde = x_plus + dx;
            
            debug_dat.cstrStateRes(n,:) = res;
            debug_dat.cstrStateKk(:,:,n) = Kk;
            
            P_tilde = P_plus;
        elseif fOpt.applyConst >= 101 && fOpt.applyConst <= 106 % fmincon
            switch (fOpt.applyConst)
                case 101
                    optimOpt.Algorithm = 'interior-point';
                    W = P_plus;
                case 102
                    optimOpt.Algorithm = 'sqp';
                    W = P_plus;
                case 103
                    optimOpt.Algorithm = 'active-set';
                    W = P_plus;
                case 104
                    optimOpt.Algorithm = 'interior-point';
                    W = I_N;
                case 105
                    optimOpt.Algorithm = 'sqp';
                    W = I_N;
                case 106
                    optimOpt.Algorithm = 'active-set';
                    W = I_N;
            end
            
            x_tilde = fmincon(@(x_tilde) L2Dist(x_tilde, x_plus, W), ...
                x_plus, [], [], [], [], [], [], ...
                @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                dLFemur, dRFemur, dLTibia, dRTibia), optimOpt);
            
            % normalize the quaternions
            x_tilde(idxOriMP) = quatnormalize(x_tilde(idxOriMP, 1)');
            x_tilde(idxOriLA) = quatnormalize(x_tilde(idxOriLA, 1)');
            x_tilde(idxOriRA) = quatnormalize(x_tilde(idxOriRA, 1)');
            P_tilde = P_plus;
        end

        xhat_con(n, :) = x_tilde;
        P_con(:, :, n)  = P_tilde;
        debug_dat.cstrState(n,:) = x_tilde;
        debug_dat.cstrP(:,:,n) = P_tilde;
        
        debug_dat.LFEO(n, :) = x_tilde(idxPosLA) + dLTibia * LTIB_CS(:, 3);
        debug_dat.RFEO(n, :) = x_tilde(idxPosRA) + dRTibia * RTIB_CS(:, 3);
        debug_dat.LFEP(n, :) = x_tilde(idxPosMP) + dPelvis/2 * PELV_CS(:, 2);
        debug_dat.RFEP(n, :) = x_tilde(idxPosMP) - dPelvis/2 * PELV_CS(:, 2);
        LFEM_z = (debug_dat.LFEP(n,:)-debug_dat.LFEO(n,:))';
        LFEM_y = LTIB_CS(:,2);
        LFEM_x = cross(LFEM_y, LFEM_z);
        RFEM_z = (debug_dat.RFEP(n,:)-debug_dat.RFEO(n,:))';
        RFEM_y = RTIB_CS(:,2);
        RFEM_x = cross(RFEM_y, RFEM_z);
        debug_dat.qLTH(n, :) = rotm2quat([LFEM_x LFEM_y LFEM_z]);
        debug_dat.qRTH(n, :) = rotm2quat([RFEM_x RFEM_y RFEM_z]);

    end
end

function y = L2Dist(x, x0, S)
%x^2 is monotomically increasing at any point not at 0, so traditional
%L2 norm involving sqrt is unnecessary, can use x^2 to find same
%location of min cost in constrained region with less computational
%cost.
%using inverse of covariance rather than I will make state est over time
%less smooth but more accurate over the average of the interval
    qW = norm(x0); %weighting for q deviation to account for diffin units between pos and q.
    res = (x-x0);
    qIdx = [7:10 17:20 27:30]';
    res(qIdx) = qW*res(qIdx);
    res = S\res; %add res*inv(S) to scale cost by certainty
    n = 1000; %have also tried 2,4,6,8,14,100,1000 around >= 14 greatly increases speed of finding solution, not much difference etween 100 and 1000
    y = (res'*res)^n;
end

function [c, ceq] = hjc_nonlcon(x, dPelvis, dLFemur, dRFemur, dLTibia, dRTibia)
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    % idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    % idxVelLA = 14:16; % column idx corresponding to the left ankle velocity
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    % idxVelRA = 24:26; % column idx corresponding to the right ankle velocity
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation

    LTIB_CS = quat2rotm(x(idxOriLA,1)');
    RTIB_CS = quat2rotm(x(idxOriRA,1)');
    PELV_CS = quat2rotm(x(idxOriMP,1)');

    LKNE = x(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA,1) + dRTibia*RTIB_CS(:,3);
    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;
    % calculate the z axis of the tibia
    LTIB_z = LTIB_CS(:,3);
    RTIB_z = RTIB_CS(:,3);

    % calculate alpha_lk and alpha_rk
    alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
    alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));

    % setup the constraint equations
    ceq = [x(idxPosLA) - x(idxPosMP) - (dPelvis/2*PELV_CS(:,2) ...
        -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
        +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
        -dLTibia*LTIB_CS(:,3)) ; ...
        x(idxPosRA) - x(idxPosMP) - (-dPelvis/2*PELV_CS(:,2)+ ...
        -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
        +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
        -dRTibia*RTIB_CS(:,3)) ];

    c = [];
end