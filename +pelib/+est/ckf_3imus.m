function [ xhat_pri, xhat_con, debug_dat ] = ckf_3imus(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, options)
% Constrained Kalman Filter for performing sensor fusion on the trajectory of
% three KMUs presumably worn on the body in the following configuration: mid
% pelvis, left ankle, right ankle
% In this state space model, the position and velocity of each kinematic
% measurement unit (KMU) is estimated in 3D space by combining the
% information from each KMU in a kalman filter. NOTE: pay special attention 
% to units:;
% position (meters)
% velocity (m/s)
% acceleration (m/2^2)
%
%   :param x0: the initial state in the GFR
%   :param gfrAccMP: the acceleration of the mid-pelvis in the GFR
%   :param gfrAccLA: the acceleration of the left ankle in the GFR
%   :param gfrAccRA: the acceleration of the right ankle in the GFR
%   :param bIsStatMP: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param bIsStatLA: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param bIsStatRA: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param qMP:		mid  pelvis orientation in the GFR (quaternion)
%   :param qLA:		left  ankle orientation in the GFR (quaternion)
%   :param qRA:		right ankle orientation in the GFR (quaternion)
%   :param dPelvis:	pelvis width
%   :param dRFemur:	right femur length
%   :param dLFemur:	left femur length
%   :param dRTibia:	right tibia length
%   :param dLTibia:	left tibia length
% 
% .. Author: - Luke Wicent Sy, Michael Del Rosario
    fOpt = struct('fs', 60, 'applyMeas', 76, 'applyCstr', 355, ...
        'sigma2QAccMP', 0.5^2, 'sigma2QAccLA', 0.5^2, 'sigma2QAccRA', 0.5^2, ...
        'sigma2RPosLA', 1e-4, 'sigma2RPosRA', 1e-4, 'sigma2RPosZMP', 1e-1, ...
        'sigma2RPosMPLimit', 1e4, 'sigma2RPosLALimit', 1e4, ...
        'sigma2RPosRALimit', 1e4, 'sigma2RPosMPLARA', 1e2, ...
        'sigma2ZuptMP', 1e-2, 'sigma2ZuptLA', 1e-2, 'sigma2ZuptRA', 1e-2, ...
        'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
        'alphaRKmin', 0, 'alphaRKmax', pi*8/9, ...
        'sckfAlpha', 0.1, 'sckfThreshold', 100, 'sckfMaxIter', 500);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    nStates = 0;
    % column idx corresponding to the mid-pelvis position
    idxPosMP = nStates+1:nStates+3; nStates = nStates + 3;
    % column idx corresponding to the mid-pelvis velocity
    idxVelMP = nStates+1:nStates+3; nStates = nStates + 3;
    % column idx corresponding to the left ankle position
    idxPosLA = nStates+1:nStates+3; nStates = nStates + 3;
    % column idx corresponding to the left ankle velocity
    idxVelLA = nStates+1:nStates+3; nStates = nStates + 3;
    % column idx corresponding to the right ankle position
    idxPosRA = nStates+1:nStates+3; nStates = nStates + 3;
    % column idx corresponding to the right ankle velocity
    idxVelRA = nStates+1:nStates+3; nStates = nStates + 3;
        
    % initialise state vector (must be column)
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', nStates});
    % demo start on floor end
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

    G = zeros(nStates,9);
    G(idxPosMP, 1:3) = dt2.*eye(3);
    G(idxVelMP, 1:3) = dt .*eye(3);
    G(idxPosLA, 4:6) = dt2.*eye(3);
    G(idxVelLA, 4:6) = dt .*eye(3);
    G(idxPosRA, 7:9) = dt2.*eye(3);
    G(idxVelRA, 7:9) = dt .*eye(3);

    % Initialise process noise covariance
    Q = diag(repelem([fOpt.sigma2QAccMP fOpt.sigma2QAccLA fOpt.sigma2QAccRA], 3));
    Q = G * Q * G';
    % initialise covariance in the state estimate
    if islogical(P0) && ~P0
        P_tilde = Q;
    elseif isscalar(P0)
        P_tilde = P0*I_N;
    else
        P_tilde = P0;
    end
       
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
    debug_dat.qLTH = nan(nSamples, 4); debug_dat.qRTH = nan(nSamples, 4);

    debug_dat.predState = nan(nSamples, nStates);
    debug_dat.predP = nan(nStates, nStates, nSamples);
    debug_dat.zuptState = nan(nSamples, nStates);
    debug_dat.zuptP = nan(nStates, nStates, nSamples);
    debug_dat.cstrState = nan(nSamples, nStates);
    debug_dat.cstrP = nan(nStates, nStates, nSamples);
    debug_dat.cstrStateU = false(nSamples, 1);
    
    debug_dat.zuptStateL = bIsStatLA;
    debug_dat.zuptStateR = bIsStatRA;
        
    nMeasure = 23;
    H = zeros(nMeasure, nStates);
    Rdiag = zeros(nMeasure, 1);
    y_k = zeros(nMeasure, nSamples);
    nMeasure = 0;
    
    % ====================================
    % applyMeas update ZUPT initialization
    % ====================================
    idxMVelMP = nMeasure+1:nMeasure+3;
    idxMVelLA = nMeasure+4:nMeasure+6;
    idxMVelRA = nMeasure+7:nMeasure+9;

    H(nMeasure+1:nMeasure+9, :) = zeros(9, nStates);
    H(idxMVelMP, idxVelMP) = eye(3);
    H(idxMVelLA, idxVelLA) = eye(3);
    H(idxMVelRA, idxVelRA) = eye(3);
    Rdiag(nMeasure+1:nMeasure+9) = repelem([fOpt.sigma2ZuptMP ...
        fOpt.sigma2ZuptLA fOpt.sigma2ZuptRA], 3);
    y_k(nMeasure+1:nMeasure+9, :) = zeros(9, nSamples);
    nMeasure = nMeasure+9;
    
    % ==============================================
    % applyMeas flat floor assumption initialization
    % ==============================================
    floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);

    idxMPosLA = nMeasure+1;     idxMPosRA = nMeasure+2;
    H(idxMPosLA, idxPosLA(3)) = eye(1);
    H(idxMPosRA, idxPosRA(3)) = eye(1);
    Rdiag(nMeasure+1:nMeasure+2) = [fOpt.sigma2RPosLA, fOpt.sigma2RPosRA];
    y_k(nMeasure+1:nMeasure+2, :) = floorZ*ones(2, nSamples);  
    nMeasure = nMeasure+2;
        
    % ===============================================
    % applyMeas covariance limiter initialization
    % ===============================================
    idxMCov1 = nMeasure+1:nMeasure+9;
    H(nMeasure+1:nMeasure+3, idxPosMP) = eye(3, 3);
    H(nMeasure+4:nMeasure+6, idxPosLA) = eye(3, 3);
    H(nMeasure+7:nMeasure+9, idxPosRA) = eye(3, 3);
    Rdiag(nMeasure+1:nMeasure+9) = repelem( [fOpt.sigma2RPosMPLimit ...
                        fOpt.sigma2RPosLALimit fOpt.sigma2RPosRALimit], 3);
    y_k(nMeasure+1:nMeasure+9, :) = zeros(9, nSamples);
    nMeasure = nMeasure+9;
    
    % ===============================================
    % applyMeas soft pelvis constraint initialization
    % ===============================================   
    % pelvis x y pos = ankle average x y pos
    idxMPosMPLARA = nMeasure+1:nMeasure+2;
    H(idxMPosMPLARA, idxPosMP(1:2)) = -eye(2,2);
    H(idxMPosMPLARA, idxPosLA(1:2)) = 0.5*eye(2,2);
    H(idxMPosMPLARA, idxPosRA(1:2)) = 0.5*eye(2,2);
    Rdiag(nMeasure+1:nMeasure+2) = repelem(fOpt.sigma2RPosMPLARA, 2);
    y_k(nMeasure+1:nMeasure+2, :) = zeros(2, nSamples);
    nMeasure = nMeasure+2;

    % pelvis z pos = initial pelvis z pos
    idxMPosZMP = nMeasure+1:nMeasure+1;
    H(idxMPosZMP, idxPosMP(3)) = 1;
    Rdiag(nMeasure+1) = fOpt.sigma2RPosZMP;
    y_k(nMeasure+1, :) = x0(idxPosMP(3));
    nMeasure = nMeasure+1;
    
    R = diag(Rdiag);
    alphalimit = struct('lkmin', fOpt.alphaLKmin, 'lkmax', fOpt.alphaLKmax, ...
                        'rkmin', fOpt.alphaRKmin, 'rkmax', fOpt.alphaRKmax);

    for n = 1:nSamples
    %% -----------------------------------------------------------------------
    % Prediction Step using accelerometer measurements ----    
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
    if fOpt.applyMeas
        idx = [];
            
        if bIsStatMP(n) idx(end+1:end+3) = idxMVelMP; end
        if bIsStatLA(n) idx(end+1:end+3) = idxMVelLA; end
        if bIsStatRA(n) idx(end+1:end+3) = idxMVelRA; end
        
        if bIsStatLA(n)
            idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
        end
        if bIsStatRA(n)
            idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
        end            
        idx(end+1:end+2) = idxMPosMPLARA;
        idx(end+1:end+1) = idxMPosZMP;
        
        res = y_k(idx, n) - H(idx, :) * x_min;
        K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx));
        x_plus = x_min + K * res;
        
        idx2 = [idx idxMCov1];
        K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        P_plus = (I_N - K * H(idx2, :)) * P_min;
        
        debug_dat.zuptState(n,:) = x_plus;
        debug_dat.zuptP(:,:,n) = P_plus;
    else
        x_plus = x_min;
        P_plus = P_min;
    end
    xhat_pos(n, :) = x_plus;
    P_pos(:, :, n)  = P_plus;
    
    %% -----------------------------------------------------------------------
    % Constraint update step ---- 
        PELV_CS = quat2rotm(qMP(n,:));
        LTIB_CS = quat2rotm(qLA(n,:));
        RTIB_CS = quat2rotm(qRA(n,:));
        % Test frankenstein constraint
        if fOpt.applyCstr
            sckfAlpha = fOpt.sckfAlpha;
            sckfThreshold = fOpt.sckfThreshold;

            x_tilde = x_plus;
            P_tilde = P_plus;

            idxCPosMP = idxPosMP;
            idxCPosLA = idxPosLA;
            idxCPosRA = idxPosRA;
            I_N2 = I_N;
            DbaseRowN = 4;
            
            % preprocessing for knee angle inequality
            % additional process to ensure knee angle does not increase
            % during iteration
            LFEM_z = x_tilde(idxCPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                     - dLTibia*LTIB_CS(:,3) - x_tilde(idxCPosLA,1);
            RFEM_z = x_tilde(idxCPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                     - dRTibia*RTIB_CS(:,3) - x_tilde(idxCPosRA,1);

            alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                             -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
            alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                             -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
            alphalimit.lkmax = min(alpha_lk, fOpt.alphaLKmax);
            alphalimit.rkmax = min(alpha_rk, fOpt.alphaRKmax);
            
            tmpLKmin = sin(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,1) - ...
                       cos(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,3);
            tmpLKmax = sin(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,1) - ...
                       cos(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,3);
            tmpRKmin = sin(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,1) - ...
                       cos(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,3);
            tmpRKmax = sin(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,1) - ...
                       cos(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,3);
                   
            for i=0:fOpt.sckfMaxIter
                % Step 1: construct D matrix              
                % calculate the z axis of femur and tibia
                LFEM_z = x_tilde(idxCPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxCPosLA,1);
                RFEM_z = x_tilde(idxCPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxCPosRA,1);
                g_dlfem = norm(LFEM_z, 2);
                g_drfem = norm(RFEM_z, 2);
                
                % D matrix construction: femur length constraint and knee hinge joint
                D = zeros(DbaseRowN, length(x_tilde));
                res = zeros(DbaseRowN, 1);
                sigmaCstr = zeros(DbaseRowN, 1);
                D(1, idxCPosMP) = LFEM_z'/g_dlfem;
                D(1, idxCPosLA) = -LFEM_z'/g_dlfem;
                D(2, idxCPosMP) = RFEM_z'/g_drfem;
                D(2, idxCPosRA) = -RFEM_z'/g_drfem;
                D(3, idxCPosMP) = LTIB_CS(:,2)';
                D(3, idxCPosLA) = -LTIB_CS(:,2)';
                D(4, idxCPosMP) = RTIB_CS(:,2)';
                D(4, idxCPosRA) = -RTIB_CS(:,2)';
                
                res(1, 1) = dLFemur - g_dlfem;
                res(2, 1) = dRFemur - g_drfem;
                res(3, 1) = (-dPelvis/2*PELV_CS(:,2) + ...
                    dLTibia*LTIB_CS(:,3))'*LTIB_CS(:,2) - D(3,:)*x_tilde;
                res(4, 1) = (dPelvis/2*PELV_CS(:,2) + ...
                    dRTibia*RTIB_CS(:,3))'*RTIB_CS(:,2) - D(4,:)*x_tilde;
                  
                % add knee inequality constraint via active set
                alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                                 -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
                alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                                 -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
                if alpha_lk < alphalimit.lkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpLKmin';
                    D(end, idxCPosLA) = -tmpLKmin';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmin ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                elseif alpha_lk > alphalimit.lkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpLKmax';
                    D(end, idxCPosLA) = -tmpLKmax';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmax ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                end
                    
                if alpha_rk < alphalimit.rkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpRKmin';
                    D(end, idxCPosRA) = -tmpRKmin';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmin ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                elseif alpha_rk > alphalimit.rkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpRKmax';
                    D(end, idxCPosRA) = -tmpRKmax';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmax ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                end
                                
                P_custom = P_tilde;
                
                % Step 3: calculate Ri               
                Ri = sckfAlpha*D*P_tilde*D'*exp(-i);
                Si_denom = diag(D * P_tilde* D');
                Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ Si_denom;
                Si_ignoreindex = (Si_denom < sigmaCstr) | isnan(Si) | (Si_denom < 1e-10);
                Si(Si_ignoreindex) = sckfThreshold+1;
                if sum(Si < sckfThreshold) == 0
                    break
                end

                % Step 4: update x_tilde and P_tilde
                Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                P_tilde = (I_N2-Kk*D)*P_tilde*(I_N2-Kk*D)' + Kk*Ri*Kk';
                dx = Kk*(res);
                x_tilde = x_tilde + dx;
            end
			
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
        LFEM_z = LFEM_z / norm(LFEM_z);
        LFEM_y = LFEM_y / norm(LFEM_y);
        LFEM_x = LFEM_x / norm(LFEM_x);
        RFEM_z = (debug_dat.RFEP(n,:)-debug_dat.RFEO(n,:))';
        RFEM_y = RTIB_CS(:,2);
        RFEM_x = cross(RFEM_y, RFEM_z);
        RFEM_y = RFEM_y / norm(RFEM_y);
        RFEM_z = RFEM_z / norm(RFEM_z);
        RFEM_x = RFEM_x / norm(RFEM_x);
        debug_dat.qLTH(n, :) = rotm2quat([LFEM_x LFEM_y LFEM_z]);
        debug_dat.qRTH(n, :) = rotm2quat([RFEM_x RFEM_y RFEM_z]);
    end
    
    debug_dat.y_k = y_k';
end