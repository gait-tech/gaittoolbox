%> @brief Constrained Lie group based extended Kalman filter implementation
%> @author Luke Wicent Sy
%> @date 14 May 2019
%>
%> EKF using Lie group/algebra representation
%> 3 KMUs presumably worn on the body in the following configuration: 
%> mid pelvis, left ankle, right ankle
%> 
%> Coding notation based on <a href="http://paulfurgale.info/news/2014/6/9/representing-robot-pose-the-good-the-bad-and-the-ugly">link</a>
%>
%> Inputs:
%>   x0       - the initial state in the GFR
%>   P0       - initial covariance
%>   bodyAccMP - the acceleration of the mid-pelvis in the body frame
%>   bodyAccLA - the acceleration of the left ankle in the body frame
%>   bodyAccRA - the acceleration of the right ankle in the body frame
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
%>   wbodyMP   - pelvis      angular velocity in the body frame
%>   wbodyLA   - left  shank angular velocity in the body frame
%>   wbodyRA   - right shank angular velocity in the body frame
%>   dPelvis   - pelvis width
%>   dRFemur   - right femur length
%>   dLFemur   - left femur length
%>   dRTibia   - right tibia length
%>   dLTibia   - left tibia length
%>   uwb_mea   - a structure containing the range measurements (m) between
%>   vel0      - struct of MP, LA, RA containing initial velocity of joints
%>   options   - struct containing the ff. settings:
%>   + fs      - sampling frequency of the magnetic and inertial measurement units
function [ xhatPri, xtilde, debug_dat ] = lieekf_3_kmus_v1(x0, P0, ...
    bodyAccMP, bIsStatMP, qMP, wbodyMP, ...
    bodyAccLA, bIsStatLA, qLA, wbodyLA, ...
    bodyAccRA, bIsStatRA, qRA, wbodyRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)
    [nSamples, ~] = size(bodyAccMP);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'sigmaQPosMP', 1, 'sigmaQPosLA', 1, 'sigmaQPosRA', 1, ...
          'sigmaQOriMP', 1, 'sigmaQOriLA', 1, 'sigmaQOriRA', 1, ...
          'sigmaQVelMP', 1, 'sigmaQVelLA', 1, 'sigmaQVelRA', 1, ...
          'sigmaQAngVelMP', 1, 'sigmaQAngVelLA', 1, 'sigmaQAngVelRA', 1, ...
          'sigmaRZPosRPV', 1e-1, 'sigmaRZPosLSK', 1e-2, 'sigmaRZPosRSK', 1e-2, ...
          'sigmaRZuptLA', 1e-1, 'sigmaRZuptRA', 1e-1);
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    addpath('liese3lib');
    
    step = struct('RPV', bIsStatMP, ...
                  'LSK', bIsStatLA, 'RSK', bIsStatRA);
    qOri = struct('RPV', qMP, 'LSK', qLA, 'RSK', qRA);
    % wOri = struct('RPV', wMP, 'LSK', wLA, 'RSK', wRA);
    u = [bodyAccMP'; bodyAccLA'; bodyAccRA';
         wbodyMP'; wbodyLA'; wbodyRA'];
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
%     xhat_pos.RPV = zeros(4,4,nSamples);
%     xhat_pos.LSK = zeros(4,4,nSamples);
%     xhat_pos.RSK = zeros(4,4,nSamples);
%     xhat_pos.RPV(1:3,1:3,:) = quat2rotm(qMP); 
%     xhat_pos.RPV(4,4,:) = 1;
%     xhat_pos.LSK(1:3,1:3,:) = quat2rotm(qLA); 
%     xhat_pos.LSK(4,4,:) = 1;
%     xhat_pos.RSK(1:3,1:3,:) = quat2rotm(qRA); 
%     xhat_pos.RSK(4,4,:) = 1;

    %% State and error covariance initialization
    bodyList = {'RPV', 'LSK', 'RSK'};
    se3StateList = {'W_T_RPV', 'W_T_LSK', 'W_T_RSK'};
    se3StateListN = length(se3StateList);
    se3StateN = se3StateListN*6;    
    r3StateN = 18;
    stateN = se3StateN + r3StateN;
    
    idx = struct('W_T_RPV', 1:6, 'W_T_LSK', 7:12, 'W_T_RSK', 13:18, ...
                 'vec', (se3StateN+1):stateN);
    r3StateList = {'vecVelMP', 'vecVelLA', 'vecVelRA', ...
                   'vecAVelMP', 'vecAVelLA', 'vecAVelRA'};
    for i=1:length(r3StateList)
        idx.(r3StateList{i}) = (1:3) + ((i-1)*3);
    end
    I_N = eye(stateN);
    
    xhatPri = struct();     % prediction update state
    xhatPos = struct();     % measurement update state
    xtilde = struct();      % constraint update state
    for i=1:se3StateListN
        sname = se3StateList{i};
        xhatPri.(sname) = nan(4,4,nSamples+1);
        xhatPos.(sname) = nan(4,4,nSamples+1);
        xtilde.(sname) = nan(4,4,nSamples+1);
    end
    xhatPri.vec = nan(r3StateN,nSamples+1);
    xhatPos.vec = nan(r3StateN,nSamples+1);
    xtilde.vec = nan(r3StateN,nSamples+1);
    
    PhatPri = nan(stateN,stateN,nSamples+1); % prediction update covariance
    PhatPos = nan(stateN,stateN,nSamples+1); % measurement update covariance
    Ptilde = nan(stateN,stateN,nSamples+1); % constraint update covariance
    
    if isscalar(P0)
        Ptilde(:,:,1) = I_N*P0;
    else
        Ptilde(:,:,1) = P0;
    end
    
    Fvec = [eye(9,9) zeros(9,9);
         zeros(9,18)];
    Gvec = [dt*eye(9,9) zeros(9,9);
         zeros(9,9) eye(9,9)];
        
    Q = diag(repelem([fOpt.sigmaQPosMP.^2, fOpt.sigmaQOriMP.^2, ...
                      fOpt.sigmaQPosLA.^2, fOpt.sigmaQOriLA.^2, ...
                      fOpt.sigmaQPosRA.^2, fOpt.sigmaQOriRA.^2, ...
                      fOpt.sigmaQVelMP.^2, fOpt.sigmaQVelLA.^2, ...
                      fOpt.sigmaQVelRA.^2, fOpt.sigmaQAngVelMP.^2, ...
                      fOpt.sigmaQAngVelLA.^2, fOpt.sigmaQAngVelRA.^2], 3));
                  
    %% H matrix initialization
    H0 = {}; y0 = {}; R0 = {};
    
    % ZUPT
    H0.lzupt = zeros(6, r3StateN);
    H0.lzupt(:,[idx.vecVelLA idx.vecAVelLA]) = eye(6,6);
    y0.lzupt = zeros(6, 1);
    R0.lzupt = repelem(fOpt.sigmaRZuptLA, 6);
    H0.rzupt = zeros(6, r3StateN);
    H0.rzupt(:,[idx.vecVelRA idx.vecAVelRA]) = eye(6,6);
    y0.rzupt = zeros(6, 1);
    R0.rzupt = repelem(fOpt.sigmaRZuptRA, 6);
    
    % zpos = some value assumption (e.g., flat floor)
    % H.zposMP, H.zposLS, H.zposRS varies with t=k
    y0.zpos = {}; R0.zpos = {};
    for i=1:se3StateListN
        sname = se3StateList{i};
        fOptparam = sprintf('sigmaRZPos%s', sname(end-2:end));
        y0.zpos.(sname) = 0;
        R0.zpos.(sname) = fOpt.(fOptparam);
    end   
    
    %% Set k=1 states (initial state taken from input)
    % SE(3) state initilization (i.e., position and orientation)
    x0Offset = [0, 10, 20];
    for i=1:se3StateListN
        sname = se3StateList{i};
        offset = x0Offset(i);
        xtilde.(sname)(:,:,1) = eye(4,4);
        xtilde.(sname)(1:3,1:3,1) = quat2rotm(x0(offset+(7:10))');
        xtilde.(sname)(1:3,4,1) = x0(offset+(1:3));
    end
    % R state initialization (i.e., velocity and angular velocity)
    xtilde.vec(:,1) = [x0([4:6 14:16 24:26]); zeros(9,1)];
    
    %% Iteration
    se32vecIdxs = {[1:3 10:12] [4:6 13:15] [7:9 16:18]};
    xiMatDebug = zeros(4,4,nSamples);
    
    for k=2:(nSamples+1)
        %% Prediction update
        F = zeros(stateN,stateN);
        for i=1:se3StateListN
            sname = se3StateList{i};
            bname = sname(end-2:end);
            
            xi = xtilde.vec(se32vecIdxs{i},k-1);
            % convert W_vel to B_vel
            xi(1:3) = quatrotate(qOri.(bname)(k-1,:), xi(1:3)')';
            bigxi = vec2tran(xi*dt);
            bufIdx = (1:6)+(i-1)*6;
            F(idx.(sname),idx.(sname)) = tranAd(bigxi);
            xhatPri.(sname)(:,:,k) = xtilde.(sname)(:,:,k-1)*bigxi;
        end
        F(idx.vec,idx.vec) = Fvec;
        xhatPri.vec(:,k) = Fvec*xtilde.vec(:,k-1) + Gvec*u(:,k-1);
        PhatPri(:,:,k) = F*Ptilde(:,:,k-1)*F' + Q;
        
        %% Measurement update
        H = {}; y = {}; R = {};
        yhat = {}; % predicted y from xhatPri
        measIdx = {}; measN = 0;
        % Flat floor assumption
        for i=1:se3StateListN
            sname = se3StateList{i};
            bname = sname(end-2:end);
            
            if step.(bname)(k) % step detected
                % Hfloor_k = blah blah
%                 Hzpos_k = zeros(1, stateN);
%                 yzpos_k = y0.zpos.(sname);
%                 Rzpos_k = R0.zpos.(sname);
%                 % TODO
%                 yhatzpos_k = y0.zpos.(sname);
                Hzpos_k = []; yzpos_k = []; yhatzpos_k = []; Rzpos_k = [];
            else
                Hzpos_k = []; yzpos_k = []; yhatzpos_k = []; Rzpos_k = [];
            end
            
            H.(sname) = Hzpos_k; R.(sname) = Rzpos_k;
            y.(sname) = yzpos_k; yhat.(sname) = yhatzpos_k;
            measIdx.(sname) = (1:size(H.(sname), 1)) + measN;
            measN = measN + size(H.(sname), 1);
        end
        
        if step.LSK(k)
            Hlzupt_k = H0.lzupt; ylzupt_k = y0.lzupt; Rlzupt_k = R0.lzupt;
        else
            Hlzupt_k = []; ylzupt_k = []; Rlzupt_k = [];
        end
        if step.RSK(k)
            Hrzupt_k = H0.rzupt; yrzupt_k = y0.rzupt; Rrzupt_k = R0.rzupt;
        else
            Hrzupt_k = []; yrzupt_k = []; Rrzupt_k = [];
        end
        H.vec = [Hlzupt_k; Hrzupt_k];
        y.vec = [ylzupt_k; yrzupt_k];
        R.vec = [Rlzupt_k; Rrzupt_k];
        yhat.vec = H.vec*xhatPri.vec(:,k);
        measIdx.(sname) = (1:size(H.vec, 1)) + measN;
        measN = measN + size(H.vec, 1);

        H.comb = zeros(measN, stateN);
        H.comb(measIdx.W_T_RPV, idx.W_T_RPV) = H.W_T_RPV;
        H.comb(measIdx.W_T_LSK, idx.W_T_LSK) = H.W_T_LSK;
        H.comb(measIdx.W_T_RSK, idx.W_T_RSK) = H.W_T_RSK;
        H.comb(measIdx.vec, idx.vec) = H.vec;
        y.comb = [y.W_T_RPV; y.W_T_LSK; y.W_T_RSK; y.vec];
        yhat.comb = [yhat.W_T_RPV; yhat.W_T_LSK; yhat.W_T_RSK; yhat.vec];
        R.comb = diag([R.W_T_RPV R.W_T_LSK R.W_T_RSK R.vec]);
        
        K = PhatPri(:,:,k)*H.comb'/(H.comb*PhatPri(:,:,k)*H.comb' + R.comb);
        PhatPos(:,:,k) = (I_N-K*H.comb)*PhatPri(:,:,k);
        measUpt = K*(y.comb-yhat.comb);
        
        for i=1:se3StateListN
            sname = se3StateList{i};
            xhatPos.(sname)(:,:,k) = xhatPri.(sname)(:,:,k)*...
                                        vec2tran(measUpt(measIdx.(sname)));
        end
        xhatPos.vec(:,k) = xhatPri.vec(:,k) + measUpt(measIdx.vec);
        
%         pRPV = pRPV + vRPV*dt + gfrAccMP(n,:)*dt2;
%         pLSK = pLSK + vLSK*dt + gfrAccLA(n,:)*dt2;
%         pRSK = pRSK + vRSK*dt + gfrAccRA(n,:)*dt2;
%         
%         vRPV = vRPV + gfrAccMP(n,:)*dt;
%         vLSK = vLSK + gfrAccLA(n,:)*dt;
%         vRSK = vRSK + gfrAccRA(n,:)*dt;
%         
%         xhat_pos.RPV(1:3,4,n) = pRPV';
%         xhat_pos.LSK(1:3,4,n) = pLSK';
%         xhat_pos.RSK(1:3,4,n) = pRSK';

        %% Constraint update
        for i=1:se3StateListN
            sname = se3StateList{i};
            xtilde.(sname)(:,:,k) = xhatPos.(sname)(:,:,k);
        end
        xtilde.vec(:,k) = xhatPos.vec(:,k);
    end
    
    for i=1:se3StateListN
        sname = se3StateList{i};
        xhatPri.(sname) = xhatPri.(sname)(:,:,2:end);
        xhatPos.(sname) = xhatPos.(sname)(:,:,2:end);
        xtilde.(sname) = xtilde.(sname)(:,:,2:end);
    end
    xhatPri.vec = xhatPri.vec(:,2:end);
    xhatPos.vec = xhatPos.vec(:,2:end);
    xtilde.vec = xtilde.vec(:,2:end);
    
    LTIBz = squeeze(xtilde.LSK(1:3,3,:))';
    RTIBz = squeeze(xtilde.RSK(1:3,3,:))';
    PELVy = squeeze(xtilde.RPV(1:3,2,:))';
    debug_dat.LFEO = squeeze(xtilde.LSK(1:3,4,:))' + dLTibia * LTIBz;
    debug_dat.RFEO = squeeze(xtilde.RSK(1:3,4,:))' + dRTibia * RTIBz;
    debug_dat.LFEP = squeeze(xtilde.RPV(1:3,4,:))' + dPelvis/2 * PELVy;
    debug_dat.RFEP = squeeze(xtilde.RPV(1:3,4,:))' - dPelvis/2 * PELVy;
    
    R_LFEM = zeros(3,3,nSamples);
    R_RFEM = zeros(3,3,nSamples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(xtilde.LSK(1:3,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(xtilde.RSK(1:3,2,:));
    RFEM_x = cross(RFEM_y, RFEM_z);
    R_RFEM(:,3,:) = RFEM_z ./ vecnorm(RFEM_z, 2, 1);
    R_RFEM(:,2,:) = RFEM_y ./ vecnorm(RFEM_y, 2, 1);
    R_RFEM(:,1,:) = RFEM_x ./ vecnorm(RFEM_x, 2, 1);
    debug_dat.qLTH = rotm2quat(R_LFEM);
    debug_dat.qRTH = rotm2quat(R_RFEM);
end