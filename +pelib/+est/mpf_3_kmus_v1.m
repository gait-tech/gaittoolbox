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
%>   uwb_mea   - a structure containing the range measurements (m) between
%>   vel0      - struct of MP, LA, RA containing initial velocity of joints
%>   options   - struct containing the ff. settings:
function [ xhat_pri, xhat_pos, debug_dat ] = pf_3_kmus_v2(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)

    %% Check that all accelerometer measurements are equal dimensions
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', nStates});
                   
    [nSamples, ~] = size(gfrAccMP);
    validateattributes(gfrAccMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(qMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    
    %% Default configurations
    fOpt = struct('fs', 60, 'applyMeas', 1, 'applyUwb', false, ...
        'applyAccBias', false, 'applyPred', 1, ...
        'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 0, 'sigmaQOriLA', 0, 'sigmaQOriRA', 0, ...
        'sigmaQOriLK', deg2rad(1), 'sigmaQOriRK', deg2rad(1), ...
        'sigmaRPosLA', 1e0, 'sigmaRPosRA', 1e0, 'sigmaRPosZMP', 1e-1, ...
        'sigmaZuptMP', 1e-1, 'sigmaZuptLA', 1e-1, 'sigmaZuptRA', 1e-1, ...
        'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
        'alphaRKmin', 0, 'alphaRKmax', pi*8/9);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% Initialization and defining the model
    idxPosMP =  1: 3; idxVelMP =  4: 6; idxPosLA =  7: 9;
    idxVelLA = 10:12; idxPosRA = 13:15; idxVelRA = 16:18;
    idxOriMP =  1: 4; idxOriLA =  5: 8; idxOriRA =  9:12;
    idxOriLK = 13; idxOriRK = 14;
    idxMMPtoLA = 1:3; idxMMPtoRA = 4:6;
    idxMVelLA = 7:9; idxMPosZLA = 10;
    idxMVelRA = 11:13; idxMPosZRA = 14;

    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt*dt;        % local variable for readability
    
    m = struct();
    m.nxn = 14;                    % Nonlinear state dimension
    m.nxl = 18;                    % Linear state dimension
    m.nx = m.nx + m.nxn;           % State dimension
    m.ny = 6+8;                    % Measurement dimension
    
    m.fn = zeros(m.nxn, 1);
    m.fl = zeros(m.nxl, 1);
    m.An = zeros(m.nxn, m.nxn);
    m.An(idxOriLK:idxOriRK,idxOriLK:idxOriRK) = eye(2,2);
    m.Al = eye(m.nxl, m.nxl);
    m.Al(idxPosMP,idxVelMP) = dt.*eye(3); % mid pelvis
    m.Al(idxPosLA,idxVelLA) = dt.*eye(3); % left ankle
    m.Al(idxPosRA,idxVelRA) = dt.*eye(3); % right ankle
    m.Gn = eye(m.nxn, m.nxn);
    m.Gl = zeros(m.nxl,9);
    m.Gl(idxPosMP, 1:3) = dt2.*eye(3);
    m.Gl(idxVelMP, 1:3) = dt .*eye(3);
    m.Gl(idxPosLA, 4:6) = dt2.*eye(3);
    m.Gl(idxVelLA, 4:6) = dt .*eye(3);
    m.Gl(idxPosRA, 7:9) = dt2.*eye(3);
    m.Gl(idxVelRA, 7:9) = dt .*eye(3);
    
    m.x0 = x0;             % Initial state
    m.P0 = P0;          % Covariance for the initial state
    % Process noise covariance
    Qn = diag([ repelem([(fOpt.sigmaQOriMP)^2 (fOpt.sigmaQOriLA)^2 ...
                         (fOpt.sigmaQOriRA)^2], 4) ...
                (fOpt.sigmaQOriLK)^2 (fOpt.sigmaQOriRK)^2 ]);
    Ql = m.Gl*diag(repelem([(fOpt.sigmaQAccMP)^2 (fOpt.sigmaQAccLA)^2 ...
                            (fOpt.sigmaQAccRA)^2], 3))*m.Gl';
    m.Q = [Ql zeros(n.nxl, n.nxn); ...
           zeros(n.nxn, n.nxl) Qn];
    m.R  =  [repelem((fOpt.sigmaZuptLA)^2, 3) fOpt.sigmaRPosLA ...
             repelem((fOpt.sigmaZuptRA)^2, 3) fOpt.sigmaRPosRA];
         
    % Define model to be used in the MPF, see eq. (18-19).
    m.h  = zeros(m.ny,1);           % Measurement model
    m.C  = zeros(m.ny,m.nxl);
    m.C(idxMMPtoLA, idxPosMP) = -eye(3, 3);
    m.C(idxMMPtoLA, idxPosLA) =  eye(3, 3);
    m.C(idxMMPtoRA, idxPosMP) = -eye(3, 3);
    m.C(idxMMPtoRA, idxPosRA) =  eye(3, 3);
    m.C(idxMVelLA,idxVelLA) = eye(3,3); m.C(idxMPosZLA,idxPosLA(3)) = 1;
    m.C(idxMVelRA,idxVelRA) = eye(3,3); m.C(idxMPosZRA,idxPosRA(3)) = 1;
    
    debug_dat = {};
    
    % local variable assignment for readability
    u_k = [gfrAccMP, gfrAccLA, gfrAccRA];
    y_k = [qMP, qLA, qRA]';
    uwbMPLARA = [uwb_mea.left_tibia_mid_pelvis,...
                 uwb_mea.mid_pelvis_right_tibia,...
                 uwb_mea.left_tibia_right_tibia];
    dBody = struct('RPV', dPelvis, 'LTH', dLFemur, 'RTH', dRFemur, ...
                   'LSK', dLTibia, 'RSK', dRTibia);
    zFloor = mean([body0.LTIO(1,3), body0.RTIO(1,3)]);
    
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);
    
    applyPred = fOpt.applyPred;
    applyMeas = fOpt.applyMeas;
    digit2ApplyMeas = mod(idivide(int32(applyMeas), 10, 'floor'), 10);
    modTenApplyMeas = mod(applyMeas, 10);
    
    switch fOpt.applyPred
        case 1
            stateTransitionFcn = @stateTransitionFcn001;
            Q = G * diag(sigmaQAccMP) * G' + diag([zeros(10, 1); deg2rad(1)*ones(8,1)]);
        case 2
            stateTransitionFcn = @stateTransitionFcn002;
            sigma = struct('accMP', fOpt.sigmaQAccMP, ...
                           'oriLK', deg2rad(0.5), 'oriRK', deg2rad(0.5));
    end
    
    if digit2ApplyMeas == 0
        measLikelihoodFcn = @measLikelihoodFcn001;
%         R = diag([1e-2*ones(14,1); (fOpt.sigmaQAccLA)^2*ones(3,1); ...
%               (fOpt.sigmaQAccRA)^2*ones(3,1)]);
        R = diag([(fOpt.sigmaROriLSK).^2*ones(3,1); ...
                  (fOpt.sigmaROriRSK).^2*ones(3,1);
                  (fOpt.sigmaZuptLA).^2*ones(3,1); (fOpt.sigmaRPosLA).^2; ...
                  (fOpt.sigmaZuptRA).^2*ones(3,1); (fOpt.sigmaRPosRA).^2; ...
                  (fOpt.sigmaRAccLA).^2*ones(3,1)
                  (fOpt.sigmaRAccRA).^2*ones(3,1)]);
    elseif digit2ApplyMeas == 1
        measLikelihoodFcn = @measLikelihoodFcn010;
        R = diag([(fOpt.sigmaRVelLA).^2*ones(3,1); ...
                  (fOpt.sigmaRVelRA).^2*ones(3,1);
                  (fOpt.sigmaRPosLA).^2; (fOpt.sigmaRPosRA).^2; ...
                  (fOpt.sigmaRVelMP).^2*ones(3,1)]);
              
        idxVelLA = 01:03; % column idx corresponding to the left ankle velocity
        idxVelRA = 04:06; % column idx corresponding to the right ankle velocity
        kfNStates = 6;
        
        F = eye(kfNStates, kfNStates);        
        G = zeros(kfNStates,6);        
        G(idxVelLA, 1:3) = dt .*eye(3);
        G(idxVelRA, 4:6) = dt .*eye(3);
        kfQ = G*diag(repelem([(fOpt.sigmaQAccLA)^2 (fOpt.sigmaQAccRA)^2], 3))*G';
        u_k = [gfrAccLA, gfrAccRA]';
        
        idxMVelLA = 1:3;
        idxMVelRA = 4:6;
        kfNMeasure = 6;
        
        H = zeros(kfNMeasure, kfNStates);
        H(idxMVelLA, idxVelLA) = eye(3);
        H(idxMVelRA, idxVelRA) = eye(3);
        Rdiag = repelem([(fOpt.sigmaZuptLA)^2 (fOpt.sigmaZuptRA)^2], 3);
        kfR = diag(Rdiag);       
        y_k = zeros(kfNMeasure, nSamples);
        
        kfI_N = eye(kfNStates,kfNStates);
        kfxPri = zeros(kfNStates,1);
        kfxPos = [vel0.LA vel0.RA]';
        kfPPri = kfI_N;
        kfPPos = 1*kfI_N;
        
        kfxListPri = nan(nSamples, kfNStates);
        kfPListPri = nan(kfNStates, kfNStates, nSamples);
        kfxListPos = nan(nSamples, kfNStates);
        kfPListPos = nan(kfNStates, kfNStates, nSamples);
    end
    
    pf = particleFilter(stateTransitionFcn, measLikelihoodFcn);
    pf.StateEstimationMethod = 'mean';
    pf.ResamplingMethod = 'systematic';
    nParticles = 10000;
    initialize(pf, nParticles, x0, P0);

    for n = 1:nSamples
        %% Predict next position. Resample particles if necessary.
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
        %% Correct position based on the given measurement to get best estimation.
        % Actual dot position is not used. Store corrected position in data array.
        if digit2ApplyMeas == 0
            step = [bIsStatLA(n, 1) bIsStatRA(n, 1)];
            if n-1 >= 1, past1Particles = xhat_pos(n-1,:)';
            else, past1Particles = x0; end
            if n-2 >= 1, past2Particles = xhat_pos(n-2,:)';
            else, past2Particles = x0; end
        
            meas = [qLA(n, :) qRA(n, :) gfrAccLA(n, :) gfrAccRA(n,:)];
            [xhat_pos(n,:), P_pos(:,:,n)] = correct(pf, meas, fs, R, applyMeas, ...
                    dBody, step, zFloor, past1Particles, past2Particles);
        elseif digit2ApplyMeas == 1
            if n-1 >= 1, past1Particles = xhat_pos(n-1,:)';
            else, past1Particles = x0; end
            
            kfxPri = F * kfxPos + G * u_k(:,n) ;
            kfPPri = F * kfPPos * F' + kfQ;
            
            nIdx = 0;
            if bIsStatLA(n, 1), nIdx = nIdx + 3; end
            if bIsStatRA(n, 1), nIdx = nIdx + 3; end
            idx = zeros(nIdx, 1); tmpIdxR = 1;
            
            if bIsStatLA(n, 1)
                idx(tmpIdxR:tmpIdxR+2) = idxMVelLA;
                tmpIdxR = tmpIdxR+3;
            end
            if bIsStatRA(n, 1)
                idx(tmpIdxR:tmpIdxR+2) = idxMVelRA;
                tmpIdxR = tmpIdxR+3;
            end
            
            res = y_k(idx, n) - H(idx, :) * kfxPri;
            K = kfPPri * H(idx, :)' /(H(idx, :) * kfPPri * H(idx,:)' + kfR(idx, idx));
            kfxPos = kfxPri + K * res;
            kfPPos = (kfI_N - K * H(idx, :)) * kfPPri;
            
            kfxListPri(n,:) = kfxPri;
            kfPListPri(:,:,n) = kfPPri;
            kfxListPos(n,:) = kfxPos;
            kfPListPos(:,:,n) = kfPPri;
        
            meas = struct();
            if bitand(modTenApplyMeas, 1)
                meas.velLA = kfxPos(idxVelLA,:)';
                meas.velRA = kfxPos(idxVelRA,:)';
            end
            if bitand(modTenApplyMeas, 2)
                meas.posZLA = zFloor; meas.posZRA = zFloor;
            end
            if bitand(modTenApplyMeas, 4)
                meas.velMP = 0;
            end
            
            [xhat_pos(n,:), P_pos(:,:,n)] = correct(pf, meas, fs, R, ...
                    applyMeas, dBody, past1Particles);
        end
    end
end

function g = mpf(y,m,u_k,step,N)
    % index list
    idxPosMP =  1: 3; idxVelMP =  4: 6; idxPosLA =  7: 9;
    idxVelLA = 10:12; idxPosRA = 13:15; idxVelRA = 16:18;
    idxOriMP =  1: 4; idxOriLA =  5: 8; idxOriRA =  9:12;
    idxOriLK = 13; idxOriRK = 14;
    idxMMPtoLA = 1:3; idxMMPtoRA = 4:6;
    idxMVelLA = 7:9; idxMPosZLA = 10;
    idxMVelRA = 11:13; idxMPosZRA = 14;
    
    Tfinal = size(y,2);         % Number of samples (measurements)
    NrPart = N;                 % Number of particles
    xf     = zeros(m.nx, Tfinal);   % Allocate room for the filtered estimates

    % (1) Initialize
    xnp = repmat(m.x0(m.nxl+1:end),1,NrPart) + ...
        chol(m.P0(m.nxl+1:end,m.nxl+1:end))*randn(m.nxn,NrPart);  % Nonlinear states
    xlp = repmat(m.x0(1:m.nxl),1,NrPart);      % Conditionally linear Gaussian states
    Pl  = m.P0(1:m.nxl,1:m.nxl);
    Pp  = repmat(Pl,[1,1,NrPart]);             % Initial covariance matrix
    xlf = zeros(size(xlp));                    % Allocate room for the filtered quantities
    Pf  = zeros(size(Pp));                     % Allocate room for the filtered quantities
    for t=1:Tfinal
        idxN = 8;
        if step(1,t), idxN = idxN+4; end % left  step detect
        if step(2,t), idxN = idxN+4; end % right step detect
        idx = 1:idxN;
        tmpiIdx = 7;
        if step(1,t) % left step detect
            idx(tmpiIdx:tmpiIdx+3) = [idxMVelLA idxMPosZLA];
            tmpiIdx = tmpiIdx + 4;
        end 
        if step(2,t) % right step detect
            idx(tmpiIdx:tmpiIdx+3) = [idxMVelRA idxMPosZRA];
            tmpiIdx = tmpiIdx + 4;
        end
        yNow = [calcLinearHJCd(qRPV, qLSK, qRSK, alphaLK, alphaRK, dBody); ...
            ];
        yhat = m.C*xlp;
        e    = repmat(yNow,1,NrPart) - yhat;
        % (2) Compute the importance weights according to eq. (25a)
        for i=1:NrPart
          M = m.C*Pp(:,:,i)*m.C' + m.R;
          q(i) = exp(-(1/2)*(e(:,i)'*inv(M)*e(:,i)));
        end
        if(sum(q)>1e-12)                % Divergence check
          q       = q/sum(q);           % Normalize the importance weights
          xf(1,t) = sum(q.*xnp,2);      % Compute estimate for the nonlinear states
          index   = sysresample(q);     % (3) Resample
          xnp     = xnp(:,index);       % Resampled nonlinear particles
          xlp     = xlp(:,index);       % Resampled linear particles
          Pp      = Pp(:,:,index);      % Resampled covariance matrices
          xlf     = xlf(:,index);       % Resampled linear particles
          Pf      = Pf(:,:,index);      % Resampled covariance matrices     
          info    = 0;
        else   % The filter has diverged
          info = 1;
          xf   = zeros(4,Tfinal);
          disp(['Weights close to zero at t=',num2str(t),' !!!']);
          return;
        end;
        % (4a) KF MU
        for i = 1:NrPart
          M         = m.C*Pp(:,:,i)*m.C' + m.R;    % Eq. (22c)
          K         = Pp(:,:,i)*m.C'*inv(M);       % Eq. (22d)
          yhat      = [(0.1*xnp(:,i).^2).*sign(xnp(:,i)); xlp(1,i) - xlp(2,i) + xlp(3,i)];
          xlf(:,i)  = xlp(:,i) + K*(yNow - yhat);  % Eq. (22a)
          Pf(:,:,i) = Pp(:,:,i) - K*M*K';          % Eq. (22b)
        end;
        xf(2:4,t) = mean(xlf,2);    % Compute estimate for the linear states
        % (4b) PF prediction according to Eq. (25b)
        xnf = xnp;
        for i = 1:NrPart
          xnp(i) = atan(xnf(i)) + xlf(1,i) + sqrt(m.An*Pf(:,:,i)*m.An' + m.Q(1,1))*randn(1);
        end;
        % (4c) KF TU
        for i = 1:NrPart
          N         = m.An*Pf(:,:,i)*m.An' + m.Q(1,1);       % Eq. (23c)
          L         = m.Al*Pf(:,:,i)*m.An'*inv(N);           % Eq. (23d)
          z         = xnp(i) - atan(xnf(i));                 % Eq. (24a)
          xlp(:,i)  = m.Al*xlf(:,i) + L*(z - m.An*xlf(:,i)); % Eq. (23a)
          Pp(:,:,i) = m.Al*Pf(:,:,i)*m.Al' + m.Q(2:end,2:end) - L*N*L'; % Eq. (23b)
        end;
    end
    g.xf   = xf;
    g.info = info;
end

function d = calcLinearHJCd(qRPV, qLSK, qRSK, alphaLK, alphaRK, dBody)
    PELV_CS = quat2rotm(qRPV);
    LTIB_CS = quat2rotm(qLSK);
    RTIB_CS = quat2rotm(qRSK);
    % setup the constraint equations               
    d = [ (dBody.RPV/2*PELV_CS(:,2) ...
             -dBody.LTH*cos(alphaLK)*LTIB_CS(:,3) ...
             +dBody.LTH*sin(alphaLK)*LTIB_CS(:,1) ...
             -dBody.LSK*LTIB_CS(:,3)) ; ...
            (-dBody.RPV/2*PELV_CS(:,2)+ ...
             -dBody.RTH*cos(alphaRK)*RTIB_CS(:,3) ...
             +dBody.RTH*sin(alphaRK)*RTIB_CS(:,1) ...
             -dBody.RSK*RTIB_CS(:,3)) ];
end