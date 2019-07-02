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
%> More detail about options
%>      fs: sampling frequency of the magnetic and inertial measurement units
%>      applyPred: 3 digit 0YX
%>          X:   1st bit use velocity to predict position
%>               2nd bit use angular velocity to predict orientation
%>          Y=0: Calculates acc and ang vel from sensor ori
%>      applyMeas: 3 digit ZYX
%>          X: Ori always on
%>               1st bit Zupt and Ankle zpos
%>               2nd bit Update angular velocity
%>          Y:   1st bit Pelvis assumption (xy=ankle average, z=initial height)
%>               2nd bit Vel cstr Y and Z
%>      applyCstr: 3 digit ZYX
%>          X:   1st bit enforce thigh length
%>               2nd bit enforce hinge knee joint
%>               3rd bit enforce knee range of motion
%>          Y=0: No P update
%>          Y=1: P update
%>
%> @param x0 initial state in the GFR
%> @param P0 initial covariance
%> @param bodyAcc acceleration of PV, LS, RS in the body frame
%> @param step boolean vector of PV, LS, RS indicating step detection
%> @param qOri PV, LS, RS orientation in the GFR (rotm)
%> @param wbody PV, LS, RS angular velocity in the body frame
%> @param body Length of PV_d (pelvis), RT_d and LT_d (r/l femur), RS_d and LS_d (r/l tibia)
%> @param uwb_mea a structure containing the range measurements (m) between
%> @param options struct containing the estimator settings:

function [ xtilde, debug_dat ] = lieekf_3_kmus_v1(x0, P0, ...
    B_a_, step, W_R_, B_w_, body, uwb_mea, options)
    N = {};
    [N.samples, ~] = size(B_a_.PV);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'applyPred', 1, 'applyMeas', 1, 'applyCstr', 1, ...
          'sigma2QAccPV', 1, 'sigma2QAccLS', 1, 'sigma2QAccRS', 1, ...
          'sigma2QGyrPV', 1, 'sigma2QGyrLS', 1, 'sigma2QGyrRS', 1, ...
          'sigma2QPosMP', 1, 'sigma2QPosLA', 1, 'sigma2QPosRA', 1, ...
          'sigma2QOriPV', 1e1, 'sigma2QOriLS', 1e1, 'sigma2QOriRS', 1e1, ...
          'sigma2QVelMP', 1, 'sigma2QVelLA', 1, 'sigma2QVelRA', 1, ...
          'sigma2QAngVelPV', 1, 'sigma2QAngVelLS', 1, 'sigma2QAngVelRS', 1, ...
          'sigma2ROriPV', 1e-3, 'sigma2ROriLS', 1e-3, 'sigma2ROriRS', 1e-3, ...
          'sigma2RAngVelPV', 1e-3, 'sigma2RAngVelLS', 1e-3, 'sigma2RAngVelRS', 1e-3, ...
          'sigma2RZPosPV', 1e-1, 'sigma2RZPosLS', 1e-4, 'sigma2RZPosRS', 1e-4, ...
          'sigma2RZuptLA', 1e-1, 'sigma2RZuptRA', 1e-1, ...
          'sigma2RVelCstrY', 1e-2, 'sigma2RVelCstrZ', 1e-2, ...
          'sigma2RXYPosPVLSRS', 1e2, ...
          'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
          'alphaRKmin', 0, 'alphaRKmax', pi*8/9 );
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% Feature on off
    knob = struct();
    knob.apModTen = mod(fOpt.applyPred, 10);
    knob.apTenDig = mod(idivide(int32(fOpt.applyPred), 10, 'floor'), 10);
    knob.amModTen = mod(fOpt.applyMeas, 10);
    knob.amTenDig = mod(idivide(int32(fOpt.applyMeas), 10, 'floor'), 10);
    knob.acModTen = mod(fOpt.applyCstr, 10);
    knob.acTenDig = mod(idivide(int32(fOpt.applyCstr), 10, 'floor'), 10);
    
    % Prediction
    % zero velocity otherwise
    knob.pred.PosWithStateVel = bitand(knob.apModTen, 1); 
    % zero angular velocity otherwise
    knob.pred.OriWithStateAngVel = bitand(knob.apModTen, 2);
    if knob.apTenDig == 0
        knob.pred.useStateinWFrameConv = false;
    else
        knob.pred.useStateinWFrameConv = true;
    end
    
    knob.meas = struct('ori', true, 'angvel', false, ...
        'zupt', false, 'xyposPVLSRS', false, ...
        'velcstrY', false, 'velcstrZ', false);   
%>          X: Ori always on
%>               1st bit Zupt and Ankle zpos
%>               2nd bit Update angular velocity
    if bitand(knob.amModTen, 1)
        knob.meas.zupt = true;
        knob.meas.zpos.LS = true;
        knob.meas.zpos.RS = true;
    else
        knob.meas.zpos.LS = false;
        knob.meas.zpos.RS = false;
        step.LS = false(N.samples, 1);
        step.RS = false(N.samples, 1);
    end
    knob.meas.angvel = bitand(knob.amModTen, 2);

%>          Y:   1st bit Pelvis assumption (xy=ankle average, z=initial height)
%>               2nd bit Vel cstr Y and Z

    if bitand(knob.amTenDig, 1)
        knob.meas.zpos.PV = true;
        step.PV = true(N.samples, 1);
        knob.meas.xyposPVLSRS = true;
    else
        knob.meas.zpos.PV = false;
        step.PV = false(N.samples, 1);
    end
    knob.meas.velcstrY = bitand(knob.amTenDig, 2);
    knob.meas.velcstrZ = bitand(knob.amTenDig, 2);
    
    % Constraint
    knob.cstr.thighlength = bitand(knob.acModTen, 1);
    knob.cstr.hingeknee = bitand(knob.acModTen, 2);
    knob.cstr.kneerom = bitand(knob.acModTen, 4);
    knob.cstr.on = knob.cstr.thighlength | knob.cstr.hingeknee | knob.cstr.kneerom;
    if knob.acTenDig == 0
        knob.cstr.Pupdate = false;
    elseif knob.acTenDig == 1
        knob.cstr.Pupdate = true;
    end
    
    addpath('liese3lib');
    
    % wOri = struct('RPV', wMP, 'LSK', wLA, 'RSK', wRA);
    g = [0 0 9.81]';        % gravity
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
    %% State and error covariance initialization
    bodyList = {'RPV', 'LSK', 'RSK'};
    se3StateList = {'W_T_PV', 'W_T_LS', 'W_T_RS'};
    N.se3StateList = length(se3StateList);
    N.se3State = N.se3StateList*6;    
    N.r3State = 18;
    
    N.state = N.se3State + N.r3State;
    
    idx = struct('W_T_PV', 1:6, 'W_T_LS', 7:12, 'W_T_RS', 13:18, ...
                 'vec', (N.se3State+1):N.state);
    r3StateList = {'vecVelPV', 'vecVelLS', 'vecVelRS', ...
                   'vecAVelPV', 'vecAVelLS', 'vecAVelRS'};
    r3StateList2 = {'velPV', 'velLS', 'velRS', ...
                    'aVelPV', 'aVelLS', 'aVelRS'};
    for i=1:length(r3StateList)
        idx.(r3StateList{i}) = (1:3) + ((i-1)*3);
        idx.(r3StateList2{i}) = (1:3) + ((i-1)*3) + N.se3State;
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
    
    F = struct(); G = struct(); Q = struct();
    F.vec = [eye(9,9) zeros(9,9);
             zeros(9,18)];
    G.vec = [dt*eye(9,9) zeros(9,9);
         zeros(9,9) eye(9,9)];
    F.curlyC = struct();
    for i=1:N.se3StateList
        sname = se3StateList{i};
        vname = sprintf('vel%s', sname(end-1:end));
        
        F.curlyC.(sname) = zeros(6, N.state);
        F.curlyC.(sname)(1:3,idx.(vname)) = dt*eye(3,3);
        if knob.pred.OriWithStateAngVel
            wname = sprintf('aVel%s', sname(end-1:end));
            F.curlyC.(sname)(4:6,idx.(wname)) = dt*eye(3,3);
        end
    end
    G.naive = zeros(N.state, 18);
    G.naive(idx.W_T_PV(1:3), 1:3) = dt2*eye(3,3);
    G.naive(idx.W_T_LS(1:3), 4:6) = dt2*eye(3,3);
    G.naive(idx.W_T_RS(1:3), 7:9) = dt2*eye(3,3);
    if knob.pred.OriWithStateAngVel
        G.naive(idx.W_T_PV(4:6), 10:12) = dt*eye(3,3);
        G.naive(idx.W_T_LS(4:6), 13:15) = dt*eye(3,3);
        G.naive(idx.W_T_RS(4:6), 16:18) = dt*eye(3,3);
    end
    G.naive(idx.vec, 1:18) = G.vec;
%     Q = diag(repelem([fOpt.sigma2QPosMP.^2, fOpt.sigma2QOriPV.^2, ...
%                       fOpt.sigma2QPosLA.^2, fOpt.sigma2QOriLS.^2, ...
%                       fOpt.sigma2QPosRA.^2, fOpt.sigma2QOriRS.^2, ...
%                       fOpt.sigma2QVelMP.^2, fOpt.sigma2QVelLA.^2, ...
%                       fOpt.sigma2QVelRA.^2, fOpt.sigma2QAngVelPV.^2, ...
%                       fOpt.sigma2QAngVelLS.^2, fOpt.sigma2QAngVelRS.^2], 3));
    Q0 = diag(repelem([fOpt.sigma2QAccPV, fOpt.sigma2QAccLS, ...
                      fOpt.sigma2QAccRS, ...
                      fOpt.sigma2QGyrPV, fOpt.sigma2QGyrLS, ...
                      fOpt.sigma2QGyrRS], 3));
    if knob.pred.OriWithStateAngVel
        Q.naive = G.naive * Q0 * G.naive';
    else
        Q.naive = G.naive * Q0 * G.naive' + ...
            diag(repelem([0, fOpt.sigma2QOriPV, 0, fOpt.sigma2QOriLS, ...
                      0, fOpt.sigma2QOriRS, 0 0 0 0 0 0], 3));
    end   
    Q.vec = G.vec*Q0*G.vec';
    Q.comb = zeros(N.state, N.state);
    Q.comb(idx.vec,idx.vec) = Q.vec;
    
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
    
    % Orientation check
    H0.ori = {};    R0.ori = {};
    for i=1:N.se3StateList
        sname = se3StateList{i};
        H0.ori.(sname) = zeros(3, N.state);
        H0.ori.(sname)(:, idx.(sname)) = [zeros(3,3) eye(3,3)];
        fOptparam = sprintf('sigma2ROri%s', sname(end-1:end));
        R0.ori.(sname) = repelem(fOpt.(fOptparam), 3);
    end
    
    % Angular velocity
    R0.angvel = repelem([fOpt.sigma2RAngVelPV, fOpt.sigma2RAngVelLS, ...
                         fOpt.sigma2RAngVelRS], 3);
    % ZUPT
    % zero velocity and angular velocity
    if false
        idx.lzupt = [idx.velLS idx.aVelLS];
        idx.rzupt = [idx.velRS idx.aVelRS];
        idx.vecLZupt = [idx.vecVelLS idx.vecAVelLS];
        idx.vecRZupt = [idx.vecVelRS idx.vecAVelRS];
    else
        idx.lzupt = idx.velLS;
        idx.rzupt = idx.velRS;
        idx.vecLZupt = idx.vecVelLS;
        idx.vecRZupt = idx.vecVelRS;
    end
    N.zupt = length(idx.lzupt);
    H0.lzupt = zeros(N.zupt, N.state);
    H0.lzupt(:,idx.lzupt) = eye(N.zupt, N.zupt);
    y0.lzupt = zeros(N.zupt, 1);
    R0.lzupt = repelem(fOpt.sigma2RZuptLA, N.zupt);
    H0.rzupt = zeros(N.zupt, N.state);
    H0.rzupt(:,idx.rzupt) = eye(N.zupt, N.zupt);
    y0.rzupt = zeros(N.zupt, 1);
    R0.rzupt = repelem(fOpt.sigma2RZuptRA, N.zupt);
        
    % zpos = some value assumption (e.g., flat floor)
    % H.zposMP, H.zposLS, H.zposRS varies with t=k
    H0.zposAnkPCircdot = point2fs([0 0 0 1]');
    y0.zpos = {}; R0.zpos = {};
    for i=1:N.se3StateList
        sname = se3StateList{i};
        fOptparam = sprintf('sigma2RZPos%s', sname(end-1:end));
        R0.zpos.(sname) = fOpt.(fOptparam);
    end
    y0.zpos.W_T_PV = xtilde.W_T_PV(3,4,1);
    y0.zpos.W_T_LS =  min([xtilde.W_T_LS(3,4,1), xtilde.W_T_RS(3,4,1)]);
    y0.zpos.W_T_RS =  y0.zpos.W_T_LS;
    
    % pelvis = ankle x y pos
    H0.p0 = [0 0 0 1]';
    H0.p0Circdot = point2fs([0 0 0 1]');
    H0.xyposDT = [eye(2,2) zeros(2,2)];
    y0.xypos = zeros(2, 1);
    R0.xypos = repelem(fOpt.sigma2RXYPosPVLSRS, 2);
    
    % velocity constraint
    y0.velcstrY = zeros(2, 1);
    y0.velcstrZ = zeros(2, 1);
    R0.velcstrY = repelem(fOpt.sigma2RVelCstrY, 2);
    R0.velcstrZ = repelem(fOpt.sigma2RVelCstrZ, 2);
    
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
    
    % knee range of motion
    d0.lkrom = 0;
    d0.rkrom = 0;
    
    %% Iteration
    se32vecIdxs = {[1:3 10:12] [4:6 13:15] [7:9 16:18]};
    u = zeros(18, N.samples);
    for k=2:(N.samples+1)
        kPast = k-1;
        %% Prediction update
        %         u(:,kPast) = [B_a_.PV(kPast,:)'; B_a_.LS(kPast,:)'; B_a_.RS(kPast,:)'; ...
%                  B_w_.PV(kPast,:)'; B_w_.LS(kPast,:)'; B_w_.RS(kPast,:)'];
        if knob.pred.useStateinWFrameConv
            u(:,kPast) = [xtilde.W_T_PV(1:3,1:3,kPast)*B_a_.PV(kPast,:)' - g; ...
                 xtilde.W_T_LS(1:3,1:3,kPast)*B_a_.LS(kPast,:)' - g; ...
                 xtilde.W_T_RS(1:3,1:3,kPast)*B_a_.RS(kPast,:)' - g; ...
                 xtilde.W_T_PV(1:3,1:3,kPast)*B_w_.PV(kPast,:)'; ...
                 xtilde.W_T_LS(1:3,1:3,kPast)*B_w_.LS(kPast,:)'; ...
                 xtilde.W_T_RS(1:3,1:3,kPast)*B_w_.RS(kPast,:)'];
        else
            u(:,kPast) = [W_R_.PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
                 W_R_.LS(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
                 W_R_.RS(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
                 W_R_.PV(:,:,kPast)*B_w_.PV(kPast,:)'; ...
                 W_R_.LS(:,:,kPast)*B_w_.LS(kPast,:)'; ...
                 W_R_.RS(:,:,kPast)*B_w_.RS(kPast,:)'];
        end
        
        F.comb = zeros(N.state,N.state);
        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);
            
            xi = xtilde.vec(se32vecIdxs{i},kPast);
            % convert W_vel to B_vel
            if knob.pred.useStateinWFrameConv
                B_R_W = xtilde.(sname)(1:3,1:3,kPast)';
            else
                B_R_W = W_R_.(bname)(:,:,kPast)';
            end
            if knob.pred.PosWithStateVel
                xi(1:3) = B_R_W*xi(1:3);
%                 xi(1:3) = B_R_W*(xi(1:3) + dt*u(se32vecIdxs{i}(1:3),kPast));
            else
                xi(1:3) = 0;
            end
            if knob.pred.OriWithStateAngVel
                xi(4:6) = B_R_W*xi(4:6);
%                 xi(4:6) = B_w_.(bname)(kPast,:)';
            else
                xi(4:6) = 0; 
            end
            bigxi = vec2tran(xi*dt);
            bigphi = vec2jac(xi*dt);
            F.comb(idx.(sname),idx.(sname)) = tranAd(bigxi);
            F.comb(idx.(sname),:) = F.comb(idx.(sname),:) + bigphi * F.curlyC.(sname);
            xhatPri.(sname)(:,:,k) = xtilde.(sname)(:,:,kPast)*bigxi;
%             Q.comb(idx.(sname),idx.(sname)) = bigphi * ...
%                         Q.vec(se32vecIdxs{i}, se32vecIdxs{i}) * bigphi';
            Q.comb(idx.(sname),idx.(sname)) = bigphi * ...
                        Q.naive(idx.(sname), idx.(sname)) * bigphi';

        end
        F.comb(idx.vec,idx.vec) = F.vec;
        xhatPri.vec(:,k) = F.vec*xtilde.vec(:,kPast) + G.vec*u(:,kPast);
        PhatPri(:,:,k) = F.comb*Ptilde(:,:,kPast)*F.comb' + Q.comb;
        
        %% Measurement update
        H = {}; deltay = {}; R = {};
        N.meas_k = 0;
        
        H.ori_k = {};   R.ori_k = {};
        H.zpos_k = {};  R.zpos_k = {};
        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);

            % Orientation
            if knob.meas.ori
                H.ori_k.(sname) = H0.ori.(sname);
                deltay.ori_k.(sname) = real(rot2vec(xhatPri.(sname)(1:3,1:3,k)' ...
                                * W_R_.(bname)(:,:,kPast)));
                R.ori_k.(sname) = R0.ori.(sname);
            else
                H.ori_k.(sname) = []; 
                deltay.ori_k.(sname) = []; 
                R.ori_k.(sname) = [];
            end
            
            N.meas_k = N.meas_k + size(H.ori_k.(sname), 1);
            
            % zPos assumption
            if step.(bname)(kPast) % step detected
                H.zpos_k.(sname) = zeros(1, N.state);
                H.zpos_k.(sname)(1, idx.(sname)) = ...
                    [0 0 1 0] * xhatPri.(sname)(:,:,k) * H0.zposAnkPCircdot;
                deltay.zpos_k.(sname) = y0.zpos.(sname) ...
                    - [0 0 1 0]*xhatPri.(sname)(:,:,k)*[0; 0; 0; 1];
                R.zpos_k.(sname) = R0.zpos.(sname);
            else
                H.zpos_k.(sname) = []; 
                deltay.zpos_k.(sname) = []; 
                R.zpos_k.(sname) = [];
            end
                       
            N.meas_k = N.meas_k + size(H.zpos_k.(sname), 1);
        end
        
        % Angular Velocity
        if knob.meas.angvel
            H.angvel = [zeros(9,18+9) eye(9,9)];
                        
            deltay.angvel = real([ W_R_.PV(:,:,kPast) * ...
                              rot2vec(xtilde.W_T_PV(1:3,1:3,kPast)' ...
                                * W_R_.PV(:,:,kPast)) ./ dt; ...
                              W_R_.LS(:,:,kPast) * ...
                              rot2vec(xtilde.W_T_LS(1:3,1:3,kPast)' ...
                                * W_R_.LS(:,:,kPast)) ./ dt; ...
                              W_R_.RS(:,:,kPast) * ...
                              rot2vec(xtilde.W_T_RS(1:3,1:3,kPast)' ...
                                * W_R_.RS(:,:,kPast)) ./ dt ]) ...
                            - xhatPri.vec(10:18,k);
            R.angvel = R0.angvel;
        else
            H.angvel = []; deltay.angvel = []; R.angvel = [];
        end
        N.meas_k = N.meas_k + size(H.angvel, 1);
        
        % XY pos
        if knob.meas.xyposPVLSRS
            H.xypos = zeros(2, N.state);
            H.xypos(:, idx.W_T_PV) = - H0.xyposDT * xhatPri.W_T_PV(:,:,k) ...
                                     * H0.p0Circdot;
            H.xypos(:, idx.W_T_LS) = H0.xyposDT * xhatPri.W_T_LS(:,:,k) ...
                                     * H0.p0Circdot / 2;
            H.xypos(:, idx.W_T_RS) = H0.xyposDT * xhatPri.W_T_RS(:,:,k) ...
                                     * H0.p0Circdot / 2;
            deltay.xypos = y0.xypos - H0.xyposDT * ...
                ( xhatPri.W_T_LS(:,:,k) * H0.p0 / 2 ...
                  + xhatPri.W_T_RS(:,:,k) * H0.p0 / 2 ...
                  - xhatPri.W_T_PV(:,:,k) * H0.p0); 
            R.xypos = R0.xypos;
        else
            H.xypos = []; deltay.xypos = []; R.xypos = [];
        end
        N.meas_k = N.meas_k + size(H.xypos, 1);
        
        % ZUPT
        if step.LS(kPast)
            H.lzupt_k = H0.lzupt; R.lzupt_k = R0.lzupt;
            deltay.lzupt_k = y0.lzupt - xhatPri.vec(idx.vecLZupt,k); 
        else
            H.lzupt_k = []; deltay.lzupt_k = []; R.lzupt_k = [];
        end
        if step.RS(kPast)
            H.rzupt_k = H0.rzupt; R.rzupt_k = R0.rzupt;
            deltay.rzupt_k = y0.rzupt - xhatPri.vec(idx.vecRZupt,k);
        else
            H.rzupt_k = []; deltay.rzupt_k = []; R.rzupt_k = [];
        end
        N.meas_k = N.meas_k + size(H.lzupt_k, 1);
        N.meas_k = N.meas_k + size(H.rzupt_k, 1);

        % Velocity Constraint
        if knob.meas.velcstrY || knob.meas.velcstrZ
            n_LT = D0.H2CT*(xhatPri.W_T_PV(:,:,k)*D0.PV_p_LH - ...
                            xhatPri.W_T_LS(:,:,k)*D0.LS_p_LK);
            n_RT = D0.H2CT*(xhatPri.W_T_PV(:,:,k)*D0.PV_p_RH - ...
                            xhatPri.W_T_RS(:,:,k)*D0.RS_p_RK);
            W_vy_PV = body.PV_d/2*xhatPri.W_T_PV(1:3, 2, k);
        end
        
        if knob.meas.velcstrY
            W_ry_LT = xhatPri.W_T_LS(1:3, 2, k);
            W_ry_RT = xhatPri.W_T_RS(1:3, 2, k);
            W_rx_LS = xhatPri.W_T_LS(1:3, 1, k);
            W_rx_RS = xhatPri.W_T_RS(1:3, 1, k);
            
            H.velcstrY = [zeros(1,18), W_ry_LT', -W_ry_LT', zeros(1,3) ...
                (hat(W_vy_PV)*W_ry_LT)', ...
                (-hat(n_LT)*W_ry_LT+W_rx_LS)', zeros(1,3);
                zeros(1,18), W_ry_RT', zeros(1,3), -W_ry_RT' ...
                (hat(-W_vy_PV)*W_ry_RT)', ...
                zeros(1,3), (-hat(n_RT)*W_ry_RT+W_rx_RS)' ];
            R.velcstrY = R0.velcstrY;
            deltay.velcstrY = -H.velcstrY(:,idx.vec)*xhatPri.vec(:,k);
        else
            H.velcstrY = []; R.velcstrY = []; deltay.velcstrY = [];
        end
        N.meas_k = N.meas_k + size(H.velcstrY, 1);

        if knob.meas.velcstrZ
            W_vz_LS = body.LS_d*xhatPri.W_T_LS(1:3, 3, k);
            W_vz_RS = body.LS_d*xhatPri.W_T_LS(1:3, 3, k);
            
            H.velcstrZ = [zeros(1,18), n_LT', -n_LT', zeros(1,3) ...
                (hat(W_vy_PV)*n_LT)', (-hat(W_vz_LS)*n_LT)', zeros(1,3);
                zeros(1,18), n_RT', zeros(1,3), -n_RT' ...
                (hat(-W_vy_PV)*n_RT)', zeros(1,3), (-hat(W_vz_RS)*n_RT)' ];
            R.velcstrZ = R0.velcstrZ;
            deltay.velcstrZ = -H.velcstrZ(:,idx.vec)*xhatPri.vec(:,k);
        else
            H.velcstrZ = []; R.velcstrZ = []; deltay.velcstrZ = [];
        end
        N.meas_k = N.meas_k + size(H.velcstrZ, 1);
        
        if N.meas_k > 0
            H.comb = [H.ori_k.W_T_PV; H.zpos_k.W_T_PV; ...
                      H.ori_k.W_T_LS; H.zpos_k.W_T_LS; ...
                      H.ori_k.W_T_RS; H.zpos_k.W_T_RS; ...
                      H.angvel; H.xypos; H.lzupt_k; H.rzupt_k; ...
                      H.velcstrY; H.velcstrZ];
            deltay.comb = [deltay.ori_k.W_T_PV; deltay.zpos_k.W_T_PV; ...
                           deltay.ori_k.W_T_LS; deltay.zpos_k.W_T_LS; ...
                           deltay.ori_k.W_T_RS; deltay.zpos_k.W_T_RS; ...
                           deltay.angvel; deltay.xypos; ...
                           deltay.lzupt_k; deltay.rzupt_k; ...
                           deltay.velcstrY; deltay.velcstrZ];
            R.comb = diag([R.ori_k.W_T_PV R.zpos_k.W_T_PV ...
                           R.ori_k.W_T_LS R.zpos_k.W_T_LS ...
                           R.ori_k.W_T_RS R.zpos_k.W_T_RS ...
                           R.angvel R.xypos R.lzupt_k R.rzupt_k ...
                           R.velcstrY R.velcstrZ]);
                
            K = PhatPri(:,:,k)*H.comb'/(H.comb*PhatPri(:,:,k)*H.comb' + R.comb);
            bigphi = eye(N.state, N.state);
            measUpt = K*(deltay.comb);

            for i=1:N.se3StateList
                sname = se3StateList{i};
                xhatPos.(sname)(:,:,k) = xhatPri.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
                % bigphi(idx.(sname), idx.(sname)) = vec2jac(measUpt(idx.(sname)));
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k) + measUpt(idx.vec);
            PhatPos(:,:,k) = bigphi*(I_N-K*H.comb)*PhatPri(:,:,k)*bigphi';
        else
            for i=1:N.se3StateList
                sname = se3StateList{i};
                xhatPos.(sname)(:,:,k) = xhatPri.(sname)(:,:,k);
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k);
            PhatPos(:,:,k) = PhatPri(:,:,k);
        end

        %% Constraint update
        if knob.cstr.on
            n_LT = D0.H2CT*(xhatPos.W_T_PV(:,:,k)*D0.PV_p_LH - ...
                                xhatPos.W_T_LS(:,:,k)*D0.LS_p_LK);
            n_RT = D0.H2CT*(xhatPos.W_T_PV(:,:,k)*D0.PV_p_RH - ...
                                xhatPos.W_T_RS(:,:,k)*D0.RS_p_RK);
            N.cstr_k = 0;
            if knob.cstr.thighlength
                N.cstr_k = N.cstr_k + 2;
            end
            if knob.cstr.hingeknee
                N.cstr_k = N.cstr_k + 2;
            end
            if knob.cstr.kneerom
                alphaLK = atan2(-dot(n_LT, xhatPos.W_T_LS(1:3,3,k)), ...
                                -dot(n_LT, xhatPos.W_T_LS(1:3,1,k))) + 0.5*pi;
                alphaRK = atan2(-dot(n_RT, xhatPos.W_T_RS(1:3,3,k)), ...
                                -dot(n_RT, xhatPos.W_T_RS(1:3,1,k))) + 0.5*pi;
                if (alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)
                    N.cstr_k = N.cstr_k + 1;
                end
                if (alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax)
                    N.cstr_k = N.cstr_k + 1;
                end
            end
            
            D = zeros(N.cstr_k, N.state);
            d = zeros(N.cstr_k, 1);
            dhat = zeros(N.cstr_k, 1);
            N.cstr_k = 0;
            
            % thigh length constraint
            if knob.cstr.thighlength
                D(N.cstr_k+1,idx.W_T_PV) = 2*n_LT'*D0.H2CT*xhatPos.W_T_PV(:,:,k)*...
                                    D0.PV_p_LH_Circdot;
                D(N.cstr_k+1,idx.W_T_LS) = -2*n_LT'*D0.H2CT*xhatPos.W_T_LS(:,:,k)*...
                                    D0.LS_p_LK_Circdot;
                D(N.cstr_k+2,idx.W_T_PV) = 2*n_RT'*D0.H2CT*xhatPos.W_T_PV(:,:,k)*...
                                    D0.PV_p_RH_Circdot;
                D(N.cstr_k+2,idx.W_T_RS) = -2*n_RT'*D0.H2CT*xhatPos.W_T_RS(:,:,k)*...
                                    D0.RS_p_RK_Circdot;
                d(N.cstr_k+(1:2),:) = [d0.ltl; d0.rtl];
                dhat(N.cstr_k+(1:2),:) = [n_LT'*n_LT; n_RT'*n_RT];
                N.cstr_k = N.cstr_k + 2;
            end
            
            % hinge knee joint constraint
            if knob.cstr.hingeknee
                W_r_LS_y = xhatPos.W_T_LS(:,:,k) * D0.p_y;
                W_r_RS_y = xhatPos.W_T_RS(:,:,k) * D0.p_y;
                D(N.cstr_k+1,idx.W_T_PV) = W_r_LS_y' * D0.H2H ...
                                  * xhatPos.W_T_PV(:,:,k) * D0.PV_p_LH_Circdot;
                D(N.cstr_k+1,idx.W_T_LS) = - W_r_LS_y' * D0.H2H ...
                                  * xhatPos.W_T_LS(:,:,k)*D0.LS_p_LK_Circdot ...
                                  + n_LT' * D0.H2CT * xhatPos.W_T_LS(:,:,k) ...
                                  * D0.p_y_Circdot;
                D(N.cstr_k+2,idx.W_T_PV) = W_r_RS_y'*D0.H2H ...
                                  * xhatPos.W_T_PV(:,:,k) * D0.PV_p_RH_Circdot;
                D(N.cstr_k+2,idx.W_T_RS) = - W_r_RS_y'*D0.H2H ...
                                  * xhatPos.W_T_RS(:,:,k)*D0.RS_p_RK_Circdot ...
                                  + n_RT'*D0.H2CT*xhatPos.W_T_RS(:,:,k) ...
                                  * D0.p_y_Circdot;
                d(N.cstr_k+(1:2),:) = [d0.lkh; d0.rkh];
                dhat(N.cstr_k+(1:2),:) = [(W_r_LS_y)'*D0.H2C*n_LT;
                               (W_r_RS_y)'*D0.H2C*n_RT];
                N.cstr_k = N.cstr_k + 2;
            end
            
            % knee range of motion
            if knob.cstr.kneerom
                if (alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)
                    alphaLK2 = min(max(alphaLK, fOpt.alphaLKmin), fOpt.alphaLKmax);
                    a = [-sin(alphaLK2-pi/2); 0; cos(alphaLK2-pi/2); 0];
                    Ta = xhatPos.W_T_LS(:,:,k) * a;

                    D(N.cstr_k+1,idx.W_T_PV) = Ta' * D0.H2H ...
                            * xhatPos.W_T_PV(:,:,k) * D0.PV_p_LH_Circdot;
                    D(N.cstr_k+1,idx.W_T_LS) = n_LT' * D0.H2CT ...
                            * xhatPos.W_T_LS(:,:,k) * point2fs(a) ...
                        - Ta' * D0.H2H * xhatPos.W_T_LS(:,:,k) ...
                            * D0.LS_p_LK_Circdot;
                    d(N.cstr_k+1,:) = d0.lkrom;
                    dhat(N.cstr_k+1,:) = n_LT' * D0.H2CT * Ta;

                    N.cstr_k = N.cstr_k + 1;
                end
                if (alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax)
                    alphaRK2 = min(max(alphaRK, fOpt.alphaRKmin), fOpt.alphaRKmax);
                    a = [-sin(alphaRK2-pi/2); 0; cos(alphaRK2-pi/2); 0];
                    Ta = xhatPos.W_T_RS(:,:,k) * a;

                    D(N.cstr_k+1,idx.W_T_PV) = Ta' * D0.H2H ...
                            * xhatPos.W_T_PV(:,:,k) * D0.PV_p_RH_Circdot;
                    D(N.cstr_k+1,idx.W_T_RS) = n_RT' * D0.H2CT ...
                            * xhatPos.W_T_RS(:,:,k) * point2fs(a) ...
                        - Ta' * D0.H2H * xhatPos.W_T_RS(:,:,k) ...
                            * D0.RS_p_RK_Circdot;
                    d(N.cstr_k+1,:) = d0.rkrom;
                    dhat(N.cstr_k+1,:) = n_RT' * D0.H2CT * Ta;

                    N.cstr_k = N.cstr_k + 1;
                end
            end
            
            K = PhatPos(:,:,k)*D'/(D*PhatPos(:,:,k)*D');
            if knob.cstr.Pupdate
                Ptilde(:,:,k) = (I_N-K*D)*PhatPos(:,:,k);
            else
                Ptilde(:,:,k) = PhatPos(:,:,k);
            end
            measUpt = K*(d-dhat);

            for i=1:N.se3StateList
                sname = se3StateList{i};
                xtilde.(sname)(:,:,k) = xhatPos.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
            end
            xtilde.vec(:,k) = xhatPos.vec(:,k) + measUpt(idx.vec);
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
    PhatPri = PhatPri(:,:,2:end);
    PhatPos = PhatPos(:,:,2:end);
    Ptilde = Ptilde(:,:,2:end);
    
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
    
    debug_dat.u = u;
    debug_dat.xhatPri = xhatPri;
    debug_dat.xhatPos = xhatPos;
    debug_dat.xtilde = xtilde;
    debug_dat.PhatPri = PhatPri;
    debug_dat.PhatPos = PhatPos;
    debug_dat.Ptilde = Ptilde;
end