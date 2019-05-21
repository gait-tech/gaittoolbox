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
%>   body.PV_d   - pelvis width
%>   body.RT_d   - right femur length
%>   body.LT_d   - left femur length
%>   body.RS_d   - right tibia length
%>   body.LS_d   - left tibia length
%>   uwb_mea   - a structure containing the range measurements (m) between
%>   vel0      - struct of MP, LA, RA containing initial velocity of joints
%>   options   - struct containing the ff. settings:
%>   + fs      - sampling frequency of the magnetic and inertial measurement units
function [ xhatPri, xtilde, debug_dat ] = lieekf_3_kmus_v1(x0, P0, ...
    bodyAccMP, bIsStatMP, qMP, wbodyMP, ...
    bodyAccLA, bIsStatLA, qLA, wbodyLA, ...
    bodyAccRA, bIsStatRA, qRA, wbodyRA, ...
    body, uwb_mea, options)
    N = {};
    [N.samples, ~] = size(bodyAccMP);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'sigmaQPosMP', 1, 'sigmaQPosLA', 1, 'sigmaQPosRA', 1, ...
          'sigmaQOriMP', 1, 'sigmaQOriLA', 1, 'sigmaQOriRA', 1, ...
          'sigmaQVelMP', 1, 'sigmaQVelLA', 1, 'sigmaQVelRA', 1, ...
          'sigmaQAngVelMP', 1, 'sigmaQAngVelLA', 1, 'sigmaQAngVelRA', 1, ...
          'sigmaRZPosPV', 1e-1, 'sigmaRZPosLS', 1e-2, 'sigmaRZPosRS', 1e-2, ...
          'sigmaRZuptLA', 1e-1, 'sigmaRZuptRA', 1e-1);
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    addpath('liese3lib');
    
    step = struct('PV', bIsStatMP, 'LS', bIsStatLA, 'RS', bIsStatRA);
    qOri = struct('PV', qMP, 'LS', qLA, 'RS', qRA);
    % wOri = struct('RPV', wMP, 'LSK', wLA, 'RSK', wRA);
    u = [bodyAccMP'; bodyAccLA'; bodyAccRA';
         wbodyMP'; wbodyLA'; wbodyRA'];
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
%     xhat_pos.RPV = zeros(4,4,N.samples);
%     xhat_pos.LSK = zeros(4,4,N.samples);
%     xhat_pos.RSK = zeros(4,4,N.samples);
%     xhat_pos.RPV(1:3,1:3,:) = quat2rotm(qMP); 
%     xhat_pos.RPV(4,4,:) = 1;
%     xhat_pos.LSK(1:3,1:3,:) = quat2rotm(qLA); 
%     xhat_pos.LSK(4,4,:) = 1;
%     xhat_pos.RSK(1:3,1:3,:) = quat2rotm(qRA); 
%     xhat_pos.RSK(4,4,:) = 1;

    %% State and error covariance initialization
    bodyList = {'RPV', 'LSK', 'RSK'};
    se3StateList = {'W_T_PV', 'W_T_LS', 'W_T_RS'};
    N.se3StateList = length(se3StateList);
    N.se3State = N.se3StateList*6;    
    N.r3State = 18;
    N.state = N.se3State + N.r3State;
    
    idx = struct('W_T_PV', 1:6, 'W_T_LS', 7:12, 'W_T_RS', 13:18, ...
                 'vec', (N.se3State+1):N.state);
    r3StateList = {'vecVelMP', 'vecVelLA', 'vecVelRA', ...
                   'vecAVelMP', 'vecAVelLA', 'vecAVelRA'};
    for i=1:length(r3StateList)
        idx.(r3StateList{i}) = (1:3) + ((i-1)*3);
    end
    I_N = eye(N.state);
    
    xhatPri = struct();     % prediction update state
    xhatPos = struct();     % measurement update state
    xtilde = struct();      % constraint update state
    for i=1:N.se3StateList
        sname = se3StateList{i};
        xhatPri.(sname) = nan(4,4,N.samples+1);
        xhatPos.(sname) = nan(4,4,N.samples+1);
        xtilde.(sname) = nan(4,4,N.samples+1);
    end
    xhatPri.vec = nan(N.r3State,N.samples+1);
    xhatPos.vec = nan(N.r3State,N.samples+1);
    xtilde.vec = nan(N.r3State,N.samples+1);
    
    PhatPri = nan(N.state,N.state,N.samples+1); % prediction update covariance
    PhatPos = nan(N.state,N.state,N.samples+1); % measurement update covariance
    Ptilde = nan(N.state,N.state,N.samples+1); % constraint update covariance
    
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
    
    %% Set k=1 states (initial state taken from input)
    % SE(3) state initilization (i.e., position and orientation)
    x0Offset = [0, 10, 20];
    for i=1:N.se3StateList
        sname = se3StateList{i};
        offset = x0Offset(i);
        xtilde.(sname)(:,:,1) = eye(4,4);
        xtilde.(sname)(1:3,1:3,1) = quat2rotm(x0(offset+(7:10))');
        xtilde.(sname)(1:3,4,1) = x0(offset+(1:3));
    end
    % R state initialization (i.e., velocity and angular velocity)
    xtilde.vec(:,1) = [x0([4:6 14:16 24:26]); zeros(9,1)];
    
    %% H matrix initialization
    H0 = {}; y0 = {}; R0 = {};
    
    % ZUPT
    % zero velocity and angular velocity
    if false
        H0.lzupt = zeros(6, N.r3State);
        H0.lzupt(:,[idx.vecVelLA idx.vecAVelLA]) = eye(6,6);
        y0.lzupt = zeros(6, 1);
        R0.lzupt = repelem(fOpt.sigmaRZuptLA, 6);
        H0.rzupt = zeros(6, N.r3State);
        H0.rzupt(:,[idx.vecVelRA idx.vecAVelRA]) = eye(6,6);
        y0.rzupt = zeros(6, 1);
        R0.rzupt = repelem(fOpt.sigmaRZuptRA, 6);
    else
        H0.lzupt = zeros(3, N.r3State);
        H0.lzupt(:,idx.vecVelLA) = eye(3,3);
        y0.lzupt = zeros(3, 1);
        R0.lzupt = repelem(fOpt.sigmaRZuptLA, 3);
        H0.rzupt = zeros(3, N.r3State);
        H0.rzupt(:,idx.vecVelRA) = eye(3,3);
        y0.rzupt = zeros(3, 1);
        R0.rzupt = repelem(fOpt.sigmaRZuptRA, 3);
    end
    
    % zpos = some value assumption (e.g., flat floor)
    % H.zposMP, H.zposLS, H.zposRS varies with t=k
    H0.zposAnkPCircdot = point2fs([0 0 0 1]');
    y0.zpos = {}; R0.zpos = {};
    for i=1:N.se3StateList
        sname = se3StateList{i};
        fOptparam = sprintf('sigmaRZPos%s', sname(end-1:end));
        y0.zpos.(sname) = xtilde.(sname)(3,4,1);
        R0.zpos.(sname) = fOpt.(fOptparam);
    end
    
    %% D matrix initialization
    D0 = {}; d0 = {};
    % thigh length constraint
    D0.H2C = [eye(3,3); zeros(1,3)]; % homogenous to cartesian
    D0.H2CT = D0.H2C';
    D0.H2H = D0.H2C*D0.H2CT;
    D0.PV_p_LH = [0 body.PV_d/2 0 1]'; 
    D0.PV_p_RH = [0 -body.PV_d/2 0 1]';
    D0.LS_p_LK = [0 0 body.LS_d 1]';    
    D0.RS_p_RK = [0 0 body.RS_d 1]';
    D0.PV_p_LH_Circdot = point2fs(D0.PV_p_LH);
    D0.PV_p_RH_Circdot = point2fs(D0.PV_p_RH);
    D0.LS_p_LK_Circdot = point2fs(D0.LS_p_LK);
    D0.RS_p_RK_Circdot = point2fs(D0.RS_p_RK);
    d0.ltl = (body.LT_d).^2;
    d0.rtl = (body.RT_d).^2;
    
    % hinge knee joint
    D0.p_y = [0 1 0 0]'; 
    D0.p_y_Circdot = point2fs(D0.p_y);
    d0.lkh = 0;
    d0.rkh = 0;
    
    %% Iteration
    se32vecIdxs = {[1:3 10:12] [4:6 13:15] [7:9 16:18]};
    
    for k=2:(N.samples+1)
        kPast = k-1;
        %% Prediction update
        F = zeros(N.state,N.state);
        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);
            
            xi = xtilde.vec(se32vecIdxs{i},kPast);
            % convert W_vel to B_vel
            xi(1:3) = quatrotate(qOri.(bname)(kPast,:), xi(1:3)')';
            bigxi = vec2tran(xi*dt);
            F(idx.(sname),idx.(sname)) = tranAd(bigxi);
            xhatPri.(sname)(:,:,k) = xtilde.(sname)(:,:,kPast)*bigxi;
        end
        F(idx.vec,idx.vec) = Fvec;
        xhatPri.vec(:,k) = Fvec*xtilde.vec(:,kPast) + Gvec*u(:,kPast);
        PhatPri(:,:,k) = F*Ptilde(:,:,kPast)*F' + Q;
        
        %% Measurement update
        H = {}; y = {}; R = {};
        yhat = {}; % predicted y from xhatPri
        measIdx = {}; N.meas_k = 0;
        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);

            % Flat floor assumption
            if step.(bname)(kPast) % step detected
                Hzpos_k = [0 0 1 0]*xhatPri.(sname)(:,:,k)*H0.zposAnkPCircdot;
                yzpos_k = y0.zpos.(sname);
                yhatzpos_k = [0 0 1 0]*xhatPri.(sname)(:,:,k)*[0; 0; 0; 1];
                Rzpos_k = R0.zpos.(sname);
            else
                Hzpos_k = []; yzpos_k = []; yhatzpos_k = []; Rzpos_k = [];
            end
            
            H.(sname) = Hzpos_k; R.(sname) = Rzpos_k;
            y.(sname) = yzpos_k; yhat.(sname) = yhatzpos_k;
            measIdx.(sname) = (1:size(H.(sname), 1)) + N.meas_k;
            N.meas_k = N.meas_k + size(H.(sname), 1);
        end
        
        if step.LS(kPast)
            Hlzupt_k = H0.lzupt; ylzupt_k = y0.lzupt; Rlzupt_k = R0.lzupt;
        else
            Hlzupt_k = []; ylzupt_k = []; Rlzupt_k = [];
        end
        if step.RS(kPast)
            Hrzupt_k = H0.rzupt; yrzupt_k = y0.rzupt; Rrzupt_k = R0.rzupt;
        else
            Hrzupt_k = []; yrzupt_k = []; Rrzupt_k = [];
        end
        H.vec = [Hlzupt_k; Hrzupt_k];
        measIdx.vec = (1:size(H.vec, 1)) + N.meas_k;
        N.meas_k = N.meas_k + size(H.vec, 1);

        if N.meas_k > 0
            y.vec = [ylzupt_k; yrzupt_k];
            R.vec = [Rlzupt_k Rrzupt_k];
            yhat.vec = H.vec*xhatPri.vec(:,k);

            H.comb = zeros(N.meas_k, N.state);
            H.comb(measIdx.W_T_PV, idx.W_T_PV) = H.W_T_PV;
            H.comb(measIdx.W_T_LS, idx.W_T_LS) = H.W_T_LS;
            H.comb(measIdx.W_T_RS, idx.W_T_RS) = H.W_T_RS;
            H.comb(measIdx.vec, idx.vec) = H.vec;
            y.comb = [y.W_T_PV; y.W_T_LS; y.W_T_RS; y.vec];
            yhat.comb = [yhat.W_T_PV; yhat.W_T_LS; yhat.W_T_RS; yhat.vec];
            R.comb = diag([R.W_T_PV R.W_T_LS R.W_T_RS R.vec]);
                
            K = PhatPri(:,:,k)*H.comb'/(H.comb*PhatPri(:,:,k)*H.comb' + R.comb);
            PhatPos(:,:,k) = (I_N-K*H.comb)*PhatPri(:,:,k);
            measUpt = K*(y.comb-yhat.comb);

            for i=1:N.se3StateList
                sname = se3StateList{i};
                xhatPos.(sname)(:,:,k) = xhatPri.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k) + measUpt(idx.vec);
        else
            for i=1:N.se3StateList
                sname = se3StateList{i};
                xhatPos.(sname)(:,:,k) = xhatPri.(sname)(:,:,k);
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k);
            PhatPos(:,:,k) = PhatPri(:,:,k);
        end

        %% Constraint update
        N.cstr_k = 4;
        if N.cstr_k > 0
            D = zeros(N.cstr_k, N.state);
            d = zeros(N.cstr_k, 1);
            dhat = zeros(N.cstr_k, 1);
            
            % thigh length constraint
            n_LT = D0.H2CT*(xhatPos.W_T_PV(:,:,k)*D0.PV_p_LH - ...
                                xhatPos.W_T_LS(:,:,k)*D0.LS_p_LK);
            n_RT = D0.H2CT*(xhatPos.W_T_PV(:,:,k)*D0.PV_p_RH - ...
                                xhatPos.W_T_RS(:,:,k)*D0.RS_p_RK);
            D(1,idx.W_T_PV) = 2*n_LT'*D0.H2CT*xhatPos.W_T_PV(:,:,k)*...
                                D0.PV_p_LH_Circdot;
            D(1,idx.W_T_LS) = -2*n_LT'*D0.H2CT*xhatPos.W_T_LS(:,:,k)*...
                                D0.LS_p_LK_Circdot;
            D(2,idx.W_T_PV) = 2*n_RT'*D0.H2CT*xhatPos.W_T_PV(:,:,k)*...
                                D0.PV_p_RH_Circdot;
            D(2,idx.W_T_RS) = -2*n_RT'*D0.H2CT*xhatPos.W_T_RS(:,:,k)*...
                                D0.RS_p_RK_Circdot;
            d(1:2,:) = [d0.ltl; d0.rtl];
            dhat(1:2,:) = [n_LT'*n_LT; n_RT'*n_RT];
            
            % hinge knee joint constraint
%             D0.p_y = [0 1 0 0]'; 
%             D0.p_y_Circdot = point2fs(D0.p_y);
            W_r_LS_y = xhatPos.W_T_LS(:,:,k)*D0.p_y;
            W_r_RS_y = xhatPos.W_T_RS(:,:,k)*D0.p_y;
            D(3,idx.W_T_PV) = W_r_LS_y'*D0.H2H*xhatPos.W_T_PV(:,:,k)*...
                                D0.PV_p_LH_Circdot;
            D(3,idx.W_T_LS) = - W_r_LS_y'*D0.H2H*xhatPos.W_T_LS(:,:,k)*...
                                D0.LS_p_LK_Circdot ...
                              + n_LT'*D0.H2CT*xhatPos.W_T_LS(:,:,k)*...
                                D0.p_y_Circdot;
            D(4,idx.W_T_PV) = W_r_RS_y'*D0.H2H*xhatPos.W_T_PV(:,:,k)*...
                                D0.PV_p_RH_Circdot;
            D(4,idx.W_T_RS) = - W_r_RS_y'*D0.H2H*xhatPos.W_T_RS(:,:,k)*...
                                D0.RS_p_RK_Circdot ...
                              + n_RT'*D0.H2CT*xhatPos.W_T_RS(:,:,k)*...
                                D0.p_y_Circdot;
            d(3:4,:) = [d0.lkh; d0.rkh];
            dhat(3:4,:) = [(W_r_LS_y)'*D0.H2C*n_LT;
                           (W_r_RS_y)'*D0.H2C*n_RT];
            
            K = PhatPos(:,:,k)*D'/(D*PhatPos(:,:,k)*D');
            % PhatPos(:,:,k) = (I_N-K*H.comb)*PhatPri(:,:,k);
            measUpt = K*(d-dhat);

            for i=1:N.se3StateList
                sname = se3StateList{i};
                xtilde.(sname)(:,:,k) = xhatPos.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
            end
            xtilde.vec(:,k) = xhatPos.vec(:,k) + measUpt(idx.vec);
            Ptilde(:,:,k) = PhatPos(:,:,k);
        else
            for i=1:N.se3StateList
                sname = se3StateList{i};
                xtilde.(sname)(:,:,k) = xhatPos.(sname)(:,:,k);
            end
            xtilde.vec(:,k) = xhatPos.vec(:,k);
            Ptilde(:,:,k) = PhatPos(:,:,k);
        end
    end
    
    for i=1:N.se3StateList
        sname = se3StateList{i};
        xhatPri.(sname) = xhatPri.(sname)(:,:,2:end);
        xhatPos.(sname) = xhatPos.(sname)(:,:,2:end);
        xtilde.(sname) = xtilde.(sname)(:,:,2:end);
    end
    xhatPri.vec = xhatPri.vec(:,2:end);
    xhatPos.vec = xhatPos.vec(:,2:end);
    xtilde.vec = xtilde.vec(:,2:end);
    
    LTIBz = squeeze(xtilde.W_T_LS(1:3,3,:))';
    RTIBz = squeeze(xtilde.W_T_RS(1:3,3,:))';
    PELVy = squeeze(xtilde.W_T_PV(1:3,2,:))';
    debug_dat.LFEO = squeeze(xtilde.W_T_LS(1:3,4,:))' + body.LS_d * LTIBz;
    debug_dat.RFEO = squeeze(xtilde.W_T_RS(1:3,4,:))' + body.RS_d * RTIBz;
    debug_dat.LFEP = squeeze(xtilde.W_T_PV(1:3,4,:))' + body.PV_d/2 * PELVy;
    debug_dat.RFEP = squeeze(xtilde.W_T_PV(1:3,4,:))' - body.PV_d/2 * PELVy;
    
    R_LFEM = zeros(3,3,N.samples);
    R_RFEM = zeros(3,3,N.samples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(xtilde.W_T_LS(1:3,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(xtilde.W_T_RS(1:3,2,:));
    RFEM_x = cross(RFEM_y, RFEM_z);
    R_RFEM(:,3,:) = RFEM_z ./ vecnorm(RFEM_z, 2, 1);
    R_RFEM(:,2,:) = RFEM_y ./ vecnorm(RFEM_y, 2, 1);
    R_RFEM(:,1,:) = RFEM_x ./ vecnorm(RFEM_x, 2, 1);
    debug_dat.qLTH = rotm2quat(R_LFEM);
    debug_dat.qRTH = rotm2quat(R_RFEM);
end