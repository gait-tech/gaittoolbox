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
function [ xhat_pri, xhat_pos, debug_dat ] = mpf_3_kmus_v2(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)

    %% Check that all accelerometer measurements are equal dimensions
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', 14+18});
                   
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
    m.nx = m.nxl + m.nxn;           % State dimension
    m.ny = 6+8;                    % Measurement dimension
    
    I_N = eye(m.nx, m.nx);
    
    m.fn = zeros(m.nxn, 1);
    m.fl = zeros(m.nxl, 1);
    m.An = zeros(m.nxn, m.nxn);
    m.An(idxOriLK:idxOriRK,idxOriLK:idxOriRK) = eye(2,2);
    m.Al = eye(m.nxl, m.nxl);
    m.Al(idxPosMP,idxVelMP) = dt.*eye(3); % mid pelvis
    m.Al(idxPosLA,idxVelLA) = dt.*eye(3); % left ankle
    m.Al(idxPosRA,idxVelRA) = dt.*eye(3); % right ankle
    m.Gn = eye(m.nxn,12); 
    m.Gl = zeros(m.nxl,9);
    m.Gl(idxPosMP, 1:3) = dt2.*eye(3);
    m.Gl(idxVelMP, 1:3) = dt .*eye(3);
    m.Gl(idxPosLA, 4:6) = dt2.*eye(3);
    m.Gl(idxVelLA, 4:6) = dt .*eye(3);
    m.Gl(idxPosRA, 7:9) = dt2.*eye(3);
    m.Gl(idxVelRA, 7:9) = dt .*eye(3);
    
    % Process noise covariance
    Qn = diag([ repelem([(fOpt.sigmaQOriMP)^2 (fOpt.sigmaQOriLA)^2 ...
                         (fOpt.sigmaQOriRA)^2], 4) ...
                (fOpt.sigmaQOriLK)^2 (fOpt.sigmaQOriRK)^2 ]);
    Ql = m.Gl*diag(repelem([(fOpt.sigmaQAccMP)^2 (fOpt.sigmaQAccLA)^2 ...
                            (fOpt.sigmaQAccRA)^2], 3))*m.Gl';
    m.Q = [Ql zeros(m.nxl, m.nxn); ...
           zeros(m.nxn, m.nxl) Qn];
    m.R  =  [repelem(0, 6) ...
             repelem((fOpt.sigmaZuptLA)^2, 3) fOpt.sigmaRPosLA ...
             repelem((fOpt.sigmaZuptRA)^2, 3) fOpt.sigmaRPosRA];
    m.R = diag(m.R);
    m.x0 = x0;             % Initial state
    % Covariance for the initial state
    if islogical(P0) && ~P0
        m.P0 = m.Q;
    elseif isscalar(P0)
        m.P0 = P0*I_N;
        idx = [idxOriMP, idxOriLA, idxOriRA]+18;
        m.P0(idx,idx) = 0;
    else
        m.P0 = P0;
    end
    
    % Define model to be used in the MPF, see eq. (18-19).
    m.h  = zeros(m.ny,1);           % Measurement model
    m.C  = zeros(m.ny,m.nxl);
    m.C(idxMMPtoLA, idxPosMP) = -eye(3, 3);
    m.C(idxMMPtoLA, idxPosLA) =  eye(3, 3);
    m.C(idxMMPtoRA, idxPosMP) = -eye(3, 3);
    m.C(idxMMPtoRA, idxPosRA) =  eye(3, 3);
    m.C(idxMVelLA,idxVelLA) = eye(3,3); m.C(idxMPosZLA,idxPosLA(3)) = 1;
    m.C(idxMVelRA,idxVelRA) = eye(3,3); m.C(idxMPosZRA,idxPosRA(3)) = 1;

    m.dBody = struct('RPV', dPelvis, 'LTH', dLFemur, 'RTH', dRFemur, ...
                     'LSK', dLTibia, 'RSK', dRTibia);
    m.zFloor = mean([x0(idxPosLA(3)) x0(idxPosRA(3))]);
    debug_dat = {};
    
    % local variable assignment for readability
    ul_k = [gfrAccMP, gfrAccLA, gfrAccRA]';
    un_k = [qMP, qLA, qRA]';
    step = [bIsStatLA(:,1) bIsStatRA(:,1)]';
%     y_k = []';
    uwbMPLARA = [uwb_mea.left_tibia_mid_pelvis,...
                 uwb_mea.mid_pelvis_right_tibia,...
                 uwb_mea.left_tibia_right_tibia];

    N = 1000;
    g = mpf(m, ul_k, un_k, step, N)
    
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
    
end

function g = mpf(m,ul_k,un_k,step,N)
    % index list
    idxPosMP =  1: 3; idxVelMP =  4: 6; idxPosLA =  7: 9;
    idxVelLA = 10:12; idxPosRA = 13:15; idxVelRA = 16:18;
    idxOriMP =  1: 4; idxOriLA =  5: 8; idxOriRA =  9:12;
    idxOriLK = 13; idxOriRK = 14;
    idxMMPtoLA = 1:3; idxMMPtoRA = 4:6;
    idxMVelLA = 7:9; idxMPosZLA = 10;
    idxMVelRA = 11:13; idxMPosZRA = 14;
    
    Tfinal = size(ul_k,2);         % Number of samples (measurements)
    NrPart = N;                 % Number of particles
    xf     = zeros(m.nx, Tfinal);   % Allocate room for the filtered estimates

    % (1) Initialize
    idxxl = 1:m.nxl;
    idxxn = m.nxl+1:m.nx;
    idxxnB = m.nxl+[idxOriLK, idxOriRK];
    
    xnp = repmat(m.x0(idxxn),1,NrPart);
    xnp([idxOriLK, idxOriRK], :) = xnp([idxOriLK, idxOriRK], :) + ...
        chol(m.P0(idxxnB,idxxnB))*randn(length(idxxnB), NrPart);  % Nonlinear states
    xlp = repmat(m.x0(idxxl),1,NrPart);      % Conditionally linear Gaussian states
    Pl  = m.P0(idxxl,idxxl);
    Pp  = repmat(Pl,[1,1,NrPart]);             % Initial covariance matrix
    xlf = zeros(size(xlp));                    % Allocate room for the filtered quantities
    Pf  = zeros(size(Pp));                     % Allocate room for the filtered quantities
    q   = zeros(1,NrPart);
    for t=1:Tfinal
        idxN = 6;
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
        yNow = [calcLinearHJCd(xlp(idxOriMP,:), ...
                               xlp(idxOriLA,:), xlp(idxOriRA,:), ...
                               xlp(idxOriLK,:), xlp(idxOriRK,:), ...
                               m.dBody); ...
                zeros(3,NrPart); repelem(m.zFloor,1,NrPart); ...
                zeros(3,NrPart); repelem(m.zFloor,1,NrPart) ];
        C2 = m.C(idx,:); R2 = m.R(idx,idx);
        Qn = m.Q(idxxn,idxxn);
        Ql = m.Q(idxxl,idxxl);
        yhat = C2*xlp;
        e    = yNow(idx)' - yhat;
        % (2) Compute the importance weights according to eq. (25a)
        for i=1:NrPart
          M = C2*Pp(:,:,i)*C2' + R2;
          q(i) = exp(-(1/2)*(e(:,i)'/M*e(:,i)));
%           q(i) = exp(-(1/2)*(e(:,i)'*inv(M)*e(:,i)));
        end
        if(sum(q)>1e-12)                % Divergence check
          q       = q/sum(q);           % Normalize the importance weights
          xf(idxxn,t) = sum(q.*xnp,2);      % Compute estimate for the nonlinear states
          index   = pelib.est.sysresample(q);     % (3) Resample
          xnp     = xnp(:,index);       % Resampled nonlinear particles
          xlp     = xlp(:,index);       % Resampled linear particles
          Pp      = Pp(:,:,index);      % Resampled covariance matrices
          xlf     = xlf(:,index);       % Resampled linear particles
          Pf      = Pf(:,:,index);      % Resampled covariance matrices     
          info    = 0;
        else   % The filter has diverged
          info = 1;
          xf   = zeros(m.n,Tfinal);
          disp(['Weights close to zero at t=',num2str(t),' !!!']);
          return;
        end
        % (4a) KF MU
        for i = 1:NrPart
          M         = C2*Pp(:,:,i)*C2' + R2;    % Eq. (22c)
          K         = Pp(:,:,i)*C2'/M;       % Eq. (22d)
          yhat      = C2*xlp(:,i);
          yNow      = [calcLinearHJCd(xlp(idxOriMP,i), ...
                               xlp(idxOriLA,i), xlp(idxOriRA,i), ...
                               xlp(idxOriLK,i), xlp(idxOriRK,i), ...
                               m.dBody); ...
                       zeros(3,1); m.zFloor; zeros(3,1); m.zFloor ];
          xlf(:,i)  = xlp(:,i) + K*(yNow(idx,1) - yhat);  % Eq. (22a)
          Pf(:,:,i) = Pp(:,:,i) - K*M*K';          % Eq. (22b)
        end
        xf(idxxl,t) = mean(xlf,2);    % Compute estimate for the linear states
        % (4b) PF prediction according to Eq. (25b)
        xnf = xnp;
        for i = 1:NrPart
          xnp(:,i) = m.An*xnf(:,i) + m.Gn*un_k(:,t) + Qn*randn(m.nxn,1);
        end
        % (4c) KF TU
        for i = 1:NrPart
          N         = m.An*Pf(:,:,i)*m.An' + Qn;             % Eq. (23c)
          L         = m.Al*Pf(:,:,i)*m.An'/N;                % Eq. (23d)
          z         = xnp(i) - x.fn;                         % Eq. (24a)
          xlp(:,i)  = m.Al*xlf(:,i) + m.Gl*ul_k(:,t) + L*(z - m.An*xlf(:,i)); % Eq. (23a)
          Pp(:,:,i) = m.Al*Pf(:,:,i)*m.Al' + Ql - L*N*L'; % Eq. (23b)
        end
    end
    g.xf   = xf;
    g.info = info;
end

function d = calcLinearHJCd(qRPV, qLSK, qRSK, alphaLK, alphaRK, dBody)
    n = length(alphaLK);
    if size(qRPV, 2) == n, qRPV = qRPV'; end
    if size(qLSK, 2) == n, qLSK = qLSK'; end
    if size(qRSK, 2) == n, qRSK = qRSK'; end
    
    PELV_CS = quat2rotm(qRPV);
    PELVy = squeeze(PELV_CS(:,2,:));
    LTIB_CS = quat2rotm(qLSK);
    LTIBx = squeeze(LTIB_CS(:,1,:));
    LTIBz = squeeze(LTIB_CS(:,3,:));
    RTIB_CS = quat2rotm(qRSK);
    RTIBx = squeeze(RTIB_CS(:,1,:));
    RTIBz = squeeze(RTIB_CS(:,3,:));
    
    % setup the constraint equations               
    d = [ (dBody.RPV/2*PELVy ...
             -dBody.LTH*cos(alphaLK).*LTIBz ...
             +dBody.LTH*sin(alphaLK).*LTIBx ...
             -dBody.LSK*LTIBz) ; ...
          (-dBody.RPV/2*PELVy+ ...
             -dBody.RTH*cos(alphaRK).*RTIBz ...
             +dBody.RTH*sin(alphaRK).*RTIBx ...
             -dBody.RSK*RTIBz) ];
end