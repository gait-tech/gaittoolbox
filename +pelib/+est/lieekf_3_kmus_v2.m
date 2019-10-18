% @brief Constrained Lie group based extended Kalman filter implementation
% @author Luke Wicent Sy
% @date 15 Aug 2019
%
% EKF using Lie group/algebra representation
% Representation: SO(3)^3 + R^27
% 3 KMUs presumably worn on the body in the following configuration: 
% mid pelvis, left ankle, right ankle
% 
% Coding notation based on <a href="http://paulfurgale.info/news/2014/6/9/representing-robot-pose-the-good-the-bad-and-the-ugly">link</a>
%
% More detail about options
%      fs: sampling frequency of the magnetic and inertial measurement units
%      applyPred: 3 digit 0YX
%          X:   1st bit predict position, velocity, angular velocity
%               2nd bit predict orientation
%          Y:   1st bit If 1 calculates acc and ang vel from state ori
%                       If 0 calculates acc and ang vel from sensor ori
%               2nd bit calculate angular velocity from orientation
%      applyMeas: 3 digit ZYX
%          X:   1st bit Zupt and ankle zpos
%               2nd bit apply angular velocity measurement
%               3rd bit apply orientation measurement
%          Y:   1st bit Pelvis xy=ankle average
%               2nd bit Pelvis z=initial height
%          Z:   1st bit covariance limiter
%      applyCstr: 3 digit ZYX
%          X:   1st bit enforce thigh length
%               2nd bit enforce hinge knee joint
%               3rd bit enforce knee range of motion
%          Z:   1st bit do P update
%
% :param x0 initial state in the GFR
% :param P0 initial covariance
% :param B_a_ acceleration of PV, LS, RS in the body frame
% :param step boolean vector of PV, LS, RS indicating step detection
% :param W_R_ PV, LS, RS orientation in the GFR (rotm)
% :param B_w_ PV, LS, RS angular velocity in the body frame
% :param body Length of PV_d (pelvis), RT_d and LT_d (r/l femur), RS_d and LS_d (r/l tibia)
% :param uwb_mea a structure containing the range measurements (m) between
% :param options struct containing the estimator settings:

function [ xtilde, debug_dat ] = lieekf_3_kmus_v2(x0, P0, ...
    B_a_, step, W_R_, B_w_, body, uwb_mea, options)
    N = {};
    [N.samples, ~] = size(B_a_.PV);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'applyPred', 1, 'applyMeas', 1, 'applyCstr', 1, ...
          'sigma2QAccPV', 1, 'sigma2QAccLS', 1, 'sigma2QAccRS', 1, ...
          'sigma2QGyrPV', 1, 'sigma2QGyrLS', 1, 'sigma2QGyrRS', 1, ...
          'sigma2QPosMP', 1, 'sigma2QPosLA', 1, 'sigma2QPosRA', 1, ...
          'sigma2QOriPV', 1e3, 'sigma2QOriLS', 1e3, 'sigma2QOriRS', 1e3, ...
          'sigma2QVelMP', 1, 'sigma2QVelLA', 1, 'sigma2QVelRA', 1, ...
          'sigma2QAngVelPV', 0, 'sigma2QAngVelLS', 0, 'sigma2QAngVelRS', 0, ...
          'sigma2ROriPV', 1e-3, 'sigma2ROriLS', 1e-3, 'sigma2ROriRS', 1e-3, ...
          'sigma2RAngVelPV', 1e-3, 'sigma2RAngVelLS', 1e-3, 'sigma2RAngVelRS', 1e-3, ...
          'sigma2RZPosPV', 1e0, 'sigma2RZPosLS', 1e-4, 'sigma2RZPosRS', 1e-4, ...
          'sigma2RZuptLA', 1e-2, 'sigma2RZuptRA', 1e-2, ...
          'sigma2RVelCstrY', 1e0, 'sigma2RVelCstrZ', 1e0, ...
          'sigma2RXYPosPVLSRS', 1e2, ...
          'sigma2RLimPos', 1e2, 'sigma2RLimOri', 1e1, ...
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
    knob.apOneDig = mod(fOpt.applyPred, 10);
    knob.apTenDig = mod(idivide(int32(fOpt.applyPred), 10, 'floor'), 10);
    knob.amOneDig = mod(fOpt.applyMeas, 10);
    knob.amTenDig = mod(idivide(int32(fOpt.applyMeas), 10, 'floor'), 10);
    knob.amHunDig = mod(idivide(int32(fOpt.applyMeas), 100, 'floor'), 10);
    knob.acOneDig = mod(fOpt.applyCstr, 10);
    knob.acTenDig = mod(idivide(int32(fOpt.applyCstr), 10, 'floor'), 10);
    knob.acHunDig = mod(idivide(int32(fOpt.applyCstr), 100, 'floor'), 10);
    
    % Prediction knobs
    % zero velocity otherwise
    knob.pred.posvelangvel = bitand(knob.apOneDig, 1);
    knob.pred.ori = bitand(knob.apOneDig, 2);
    knob.pred.useStateinWFrameConv = bitand(knob.apTenDig, 1);
    knob.pred.useAngvelFromOri = bitand(knob.apTenDig, 2);
    knob.pred.indepPosPred = bitand(knob.apTenDig, 4);
    
    % Measurement knobs
    knob.meas = struct('angvel', bitand(knob.amOneDig, 2), ...
        'ori', bitand(knob.amOneDig, 4), ...
        'zpos', struct('LS', false, 'RS', false, 'PV', false), ...
        'zupt', false, 'xyposPVLSRS', false, ...
        'velcstrY', false, 'velcstrZ', false, ...
        'covlim', bitand(knob.amHunDig, 1) );   

    if bitand(knob.amOneDig, 1)
        knob.meas.zupt = true;
        knob.meas.zpos.LS = true;
        knob.meas.zpos.RS = true;
    else
        step.LS = false(N.samples, 1);
        step.RS = false(N.samples, 1);
    end
    knob.meas.xyposPVLSRS = bitand(knob.amTenDig, 1);
    if bitand(knob.amTenDig, 2)
        knob.meas.zpos.PV = true;
        step.PV = true(N.samples, 1);
    else
        knob.meas.zpos.PV = false;
        step.PV = false(N.samples, 1);
    end
    
    knob.meas.velcstrYAlways = bitshift(knob.amHunDig, -1) == 1;
    knob.meas.velcstrYStep = bitshift(knob.amHunDig, -1) == 2;
%     knob.meas.velcstrYAlways = 0;
%     knob.meas.velcstrYStep = 0;
    knob.meas.velcstrZAlways = bitshift(knob.amHunDig, -1) == 1;
    knob.meas.velcstrZStep = bitshift(knob.amHunDig, -1) == 2;
%     knob.meas.velcstrZAlways = 0;
%     knob.meas.velcstrZStep = 0;
    
    % Constraint
    knob.cstr.thighlength = bitand(knob.acOneDig, 1);
    knob.cstr.hingeknee = bitand(knob.acOneDig, 2);
    knob.cstr.kneerom = bitand(knob.acOneDig, 4);
    knob.cstr.Pupdate = bitand(knob.acHunDig, 1);

    knob.cstr.velcstrYAlways = bitand(knob.acTenDig, 3) == 1;
    knob.cstr.velcstrYStep = bitand(knob.acTenDig, 3) == 2;
%     knob.cstr.velcstrYAlways = 0;
%     knob.cstr.velcstrYStep = 0;
    knob.cstr.velcstrZAlways = bitand(knob.acTenDig, 3) == 1;
    knob.cstr.velcstrZStep = bitand(knob.acTenDig, 3) == 2;
%     knob.cstr.velcstrZAlways = 0;
%     knob.cstr.velcstrZStep = 0;
    
    addpath('liese3lib');
       
    %% State and error covariance initialization
%     so3StateList = ["W_R_PV", "W_R_LS", "W_R_RS"];
%     bodyList = ["PV", "LS", "RS"];
    so3StateList = ["W_R_PV", "W_R_LS", "W_R_RS"];
    bodyList = ["PV", "LS", "RS"];
    N.so3StateList = length(so3StateList);
    N.so3State = N.so3StateList*3;    
    N.r3State = 9*3;
    N.state = N.so3State;
    
    % calculate indices
    idx = struct('so3State', 1:9, ...
                 'W_R_PV', 1:3, 'W_R_LS', 4:6, 'W_R_RS', 7:9, ...
                 'vec', struct() );
    for i=["pos", "vel", "avel"]
        for j=bodyList
            k = sprintf('%s%s', i, j);
            idx.(k) = (1:3) + N.state;
            idx.vec.(k) = (1:3) + N.state - N.so3State;
            N.state = N.state + 3;
        end
    end
    idx.vecState = (N.so3State+1):N.state;
    
    g = [0 0 9.81]';        % gravity
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    I_N = eye(N.state);
    I_3 = eye(3);
    
    xhatPri = struct();     % prediction update state
    xhatPos = struct();     % measurement update state
    xtilde = struct();      % constraint update state
    for i=so3StateList
        xhatPri.(i) = nan(3,3,N.samples+1);
        xhatPos.(i) = nan(3,3,N.samples+1);
        xtilde.(i) = nan(3,3,N.samples+1);
    end
    xhatPri.vec = nan(N.r3State,N.samples+1);
    xhatPos.vec = nan(N.r3State,N.samples+1);
    xtilde.vec = nan(N.r3State,N.samples+1);
    
    PhatPri = nan(N.state,N.state,N.samples+1); % prediction update covariance
    PhatPos = nan(N.state,N.state,N.samples+1); % measurement update covariance
    Ptilde = nan(N.state,N.state,N.samples+1); % constraint update covariance
    
    %% Set initial states (k=1)
    xtilde.W_R_PV(:,:,1) = quat2rotm(x0(07:10)');
    xtilde.W_R_LS(:,:,1) = quat2rotm(x0(17:20)');
    xtilde.W_R_RS(:,:,1) = quat2rotm(x0(27:30)');
    xtilde.vec(:,1) = [x0([01:03 11:13 21:23]); ... % position
                       x0([04:06 14:16 24:26]); ... % velocity
                       zeros(9,1)]; % angular velocity
    if isscalar(P0)
        P0 = I_N*P0;
    end
    Ptilde(:,:,1) = P0;
    
    %% Prediction step initialization
    % prediction model
    F = struct(); G = struct(); Q = struct();
    F.vec = [eye(9,9) dt*eye(9,9) zeros(9,9);
             zeros(9,9) eye(9,9) zeros(9,9);
             zeros(9,27)];
    G.vec = [dt2*eye(9,9) zeros(9,9)
             dt*eye(9,9) zeros(9,9);
             zeros(9,9) eye(9,9)];
    F.curlyC = zeros(N.state, N.state);
    F.curlyC(idx.posPV, idx.velPV) = dt*eye(3,3);
    F.curlyC(idx.posLS, idx.velLS) = dt*eye(3,3);
    F.curlyC(idx.posRS, idx.velRS) = dt*eye(3,3);
    
    Q.raw0 = diag(repelem([fOpt.sigma2QAccPV, fOpt.sigma2QAccLS, ...
                           fOpt.sigma2QAccRS, ...
                           fOpt.sigma2QGyrPV, fOpt.sigma2QGyrLS, ...
                           fOpt.sigma2QGyrRS], 3));
    Q.vec = G.vec*Q.raw0*G.vec';
    Q.comb = zeros(N.state,N.state);
    if knob.pred.ori
        Q.comb(idx.so3State,idx.so3State) = diag(repelem( ...
                [fOpt.sigma2QAngVelPV*(dt.^2), ...
                 fOpt.sigma2QAngVelLS*(dt.^2), ...
                 fOpt.sigma2QAngVelRS*(dt.^2)], 3));
    else
        Q.comb(idx.so3State,idx.so3State) = diag(repelem( ...
                [fOpt.sigma2QOriPV, fOpt.sigma2QOriLS, fOpt.sigma2QOriRS], 3));
    end
    Q.comb(idx.vecState,idx.vecState) = Q.vec;

    if knob.pred.useAngvelFromOri
        B_w2_ = struct();
        for i=bodyList
            B_w2_.(i) = zeros(N.samples, 3);
            for j=1:(N.samples-1)
                B_w2_.(i)(j, :) = rot2vec(W_R_.(i)(:,:,j)'*W_R_.(i)(:,:,j+1))/dt;
            end
        end
    else
        B_w2_ = B_w_;
    end
    
    idx.u = struct('accPV', 1:3, 'accLS', 4:6, 'accRS', 7:9, ...
                   'avelPV', 10:12, 'avelLS', 13:15, 'avelRS', 16:18);
    u = zeros(18, N.samples);
    
    %% Measurement step initialization
    H = {}; y = {}; R = {};
    
    % update orientation measurement
    H.ori = [eye(9,9) zeros(9,27)];
    R.ori = repelem([fOpt.sigma2ROriPV, fOpt.sigma2ROriLS, ...
                     fOpt.sigma2ROriRS], 3);
    % y.ori is variable and will be generated on the spot
    
    % Angular velocity
    H.angvel = [zeros(9,27) eye(9,9)];
    R.angvel = repelem([fOpt.sigma2RAngVelPV, fOpt.sigma2RAngVelLS, ...
                        fOpt.sigma2RAngVelRS], 3);
    % y.angvel is variable and will be generated on the spot
    
    % Zero velocity and angular velocity update
%     idx.lzupt = [idx.velLS idx.avelLS];
%     idx.rzupt = [idx.velRS idx.avelRS];
%     idx.vec.lzupt = [idx.vec.velLS idx.vec.avelLS];
%     idx.vec.rzupt = [idx.vec.velRS idx.vec.avelRS];
    idx.lzupt = idx.velLS;
    idx.rzupt = idx.velRS;    
    N.zupt = length(idx.lzupt);
    
    H.lzupt = zeros(N.zupt, N.state);
    H.lzupt(:, idx.lzupt) = eye(N.zupt, N.zupt);
    y.lzupt = zeros(N.zupt, 1);
    R.lzupt = repelem(fOpt.sigma2RZuptLA, N.zupt);
    
    H.rzupt = zeros(N.zupt, N.state);
    H.rzupt(:, idx.rzupt) = eye(N.zupt, N.zupt);
    y.rzupt = zeros(N.zupt, 1);
    R.rzupt = repelem(fOpt.sigma2RZuptRA, N.zupt);
        
    % zpos = some value assumption (e.g., flat floor)
    H.zpos = {}; y.zpos = {}; R.zpos = {};
    for i = bodyList
        H.zpos.(i) = zeros(1, N.state);
        j = sprintf("pos%s", i);
        H.zpos.(i)(idx.(j)(3)) = 1;
        j = sprintf("sigma2RZPos%s", i);
        R.zpos.(i) = fOpt.(j);
    end
    y.zpos.PV = xtilde.vec(idx.vec.posPV(3),1);
    y.zpos.LS = min([xtilde.vec(idx.vec.posLS(3),1), ...
                     xtilde.vec(idx.vec.posRS(3),1)]);
    y.zpos.RS = y.zpos.LS;
    
    % pelvis = ankle x y pos
    H.xypos = [zeros(2,9) ...
               [1 0 0 -0.5 0 0 -0.5 0 0;
                0 1 0 0 -0.5 0 0 -0.5 0] ...
               zeros(2,18)];
    y.xypos = zeros(2, 1);
    R.xypos = repelem(fOpt.sigma2RXYPosPVLSRS, 2);
    
    % covariance limiter
    H.covlim = [eye(18) zeros(18, N.state-18)];
    R.covlim = repelem([fOpt.sigma2RLimOri, fOpt.sigma2RLimPos], 9);
    
    %% Constraint step initialization
    D = struct('i_x', [1 0 0]', 'i_y', [0 1 0]', 'i_z', [0 0 1]');
    D.ihat_x = hat(D.i_x); D.ihat_y = hat(D.i_y); D.ihat_z = hat(D.i_z); 
    d = {};
    
    % thigh length constraint
    D.PV_p_LH = [0 body.PV_d/2 0]';     D.PV_phat_LH = hat(D.PV_p_LH);
    D.PV_p_RH = [0 -body.PV_d/2 0]';    D.PV_phat_RH = hat(D.PV_p_RH);
    D.LS_p_LK = [0 0 body.LS_d]';       D.LS_phat_LK = hat(D.LS_p_LK);  
    D.RS_p_RK = [0 0 body.RS_d]';       D.RS_phat_RK = hat(D.RS_p_RK);
    d.ltl = (body.LT_d).^2;
    d.rtl = (body.RT_d).^2;
    
    % hinge knee joint
    d.lkh = 0;
    d.rkh = 0;
    
    % knee range of motion
    d.lkrom = 0;
    d.rkrom = 0;
    
    %% Debug
    debug_dat.measUptPos = zeros(N.state, N.samples);
    debug_dat.measUptTilde = zeros(N.state, N.samples);
    
    %% Iteration
    for k=2:(N.samples+1)
        kPast = k-1;
        
        %% Prediction update
        if knob.pred.useStateinWFrameConv
            u(:,kPast) = [xtilde.W_R_PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
                 xtilde.W_R_LS(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
                 xtilde.W_R_RS(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
                 xtilde.W_R_PV(:,:,kPast)*B_w2_.PV(kPast,:)'; ...
                 xtilde.W_R_LS(:,:,kPast)*B_w2_.LS(kPast,:)'; ...
                 xtilde.W_R_RS(:,:,kPast)*B_w2_.RS(kPast,:)'];
        else
            u(:,kPast) = [W_R_.PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
                 W_R_.LS(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
                 W_R_.RS(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
                 W_R_.PV(:,:,kPast)*B_w2_.PV(kPast,:)'; ...
                 W_R_.LS(:,:,kPast)*B_w2_.LS(kPast,:)'; ...
                 W_R_.RS(:,:,kPast)*B_w2_.RS(kPast,:)'];
        end
        
        F.AdG = eye(N.state, N.state);
        bigphi = eye(N.state, N.state);
        for i=so3StateList
            j = string(i{1}(end-1:end));
            
            if knob.pred.ori
                xi = B_w2_.(j)(kPast,:)'; % B_avel
            else
                xi = zeros(3,1);
            end
            
            xhatPri.(i)(:,:,k) = xtilde.(i)(:,:,kPast)*vec2rot(xi*dt);
            
            bigphi(idx.(i), idx.(i)) = vec2jac(-xi*dt);
            % for SO(3) Ad(X) = X;
            % F.AdG(idx.(i), idx.(i)) = tranAd(vec2rot(-xi*dt));
            F.AdG(idx.(i), idx.(i)) = vec2rot(-xi*dt);
        end
        
        if knob.pred.posvelangvel
            xhatPri.vec(:,k) = F.vec*xtilde.vec(:,kPast) + G.vec*u(:,kPast);
        else
            xhatPri.vec(:,k) = xtilde.vec(:,kPast);
        end
        
        F.comb = F.AdG + bigphi*F.curlyC;
        PhatPri(:,:,k) = F.comb*Ptilde(:,:,kPast)*F.comb' + bigphi*Q.comb*bigphi';
        
        %% Measurement update
        N.meas_k = (knob.meas.ori > 0) * 9 + ... % orientation
                   (knob.meas.zpos.PV && step.PV(kPast)) + ... %zpos PV
                   (knob.meas.zpos.LS && step.LS(kPast)) + ... %zpos LS
                   (knob.meas.zpos.RS && step.RS(kPast)) + ... %zpos RS
                   (knob.meas.angvel > 0) * 9 + ... % angular velocity
                   (knob.meas.xyposPVLSRS > 0) * 2 + ... % PV = ankle xy pos
                   (knob.meas.zupt && step.LS(kPast)) * 3 + ... %lzupt
                   (knob.meas.zupt && step.RS(kPast)) * 3 ; %rzupt
          
        H.comb = zeros(N.meas_k, N.state);
        dy = zeros(N.meas_k, 1);
        R.comb = zeros(N.meas_k, 1);
        N.meas_k = 0;

        % Update orientation from measurement
        if knob.meas.ori
            H.comb(N.meas_k+1:N.meas_k+9, :) = H.ori;
            R.comb(N.meas_k+1:N.meas_k+9, :) = R.ori;
            dy(N.meas_k+1:N.meas_k+9, :) = [ ...
                real(rot2vec(xhatPri.W_R_PV(:,:,k)' * W_R_.PV(:,:,kPast))); ...
                real(rot2vec(xhatPri.W_R_LS(:,:,k)' * W_R_.LS(:,:,kPast))); ...
                real(rot2vec(xhatPri.W_R_RS(:,:,k)' * W_R_.RS(:,:,kPast))) ];
            N.meas_k = N.meas_k + 9;
        end
        
        % z pos assumption
        for i=bodyList
            if knob.meas.zpos.(i) && step.(i)(kPast) % step detected
                H.comb(N.meas_k+1, :) = H.zpos.(i);
                R.comb(N.meas_k+1, :) = R.zpos.(i);
                j = sprintf('pos%s', i);
                dy(N.meas_k+1, :) = y.zpos.(i) - xhatPri.vec(idx.vec.(j)(3),k);
                N.meas_k = N.meas_k + 1;
            end
        end
        
        % Update angular velocity from measurement
        if knob.meas.angvel
            H.comb(N.meas_k+1:N.meas_k+9, :) = H.angvel;
            R.comb(N.meas_k+1:N.meas_k+9, :) = R.angvel;
            dy(N.meas_k+1:N.meas_k+9, :) = [ ...
                u(idx.u.avelPV,kPast) - xhatPri.vec(idx.vec.avelPV,k); ...
                u(idx.u.avelLS,kPast) - xhatPri.vec(idx.vec.avelLS,k); ...
                u(idx.u.avelRS,kPast) - xhatPri.vec(idx.vec.avelRS,k) ];
            N.meas_k = N.meas_k + 9;
        end
        
        % PV = ankle XY position assumption
        if knob.meas.xyposPVLSRS
            H.comb(N.meas_k+1:N.meas_k+2, :) = H.xypos;
            R.comb(N.meas_k+1:N.meas_k+2, :) = R.xypos;
            dy(N.meas_k+1:N.meas_k+2, :) = y.xypos - ...
                H.xypos(:,idx.vecState)*xhatPri.vec(:,k);
            N.meas_k = N.meas_k + 2;
        end
        
        % ZUPT
        if knob.meas.zupt && step.LS(kPast)
            H.comb(N.meas_k+1:N.meas_k+3, :) = H.lzupt;
            R.comb(N.meas_k+1:N.meas_k+3, :) = R.lzupt;
            dy(N.meas_k+1:N.meas_k+3, :) = y.lzupt - xhatPri.vec(idx.vec.velLS,k);
            N.meas_k = N.meas_k + 3;
        end
        if knob.meas.zupt && step.RS(kPast)
            H.comb(N.meas_k+1:N.meas_k+3, :) = H.rzupt;
            R.comb(N.meas_k+1:N.meas_k+3, :) = R.rzupt;
            dy(N.meas_k+1:N.meas_k+3, :) = y.rzupt - xhatPri.vec(idx.vec.velRS,k);
            N.meas_k = N.meas_k + 3;
        end
        
        % apply measurement update
        if N.meas_k > 0               
            K = PhatPri(:,:,k)*H.comb'/(H.comb*PhatPri(:,:,k)*H.comb' + diag(R.comb));
            bigphi = eye(N.state, N.state);
            measUpt = K*dy;
            debug_dat.measUptPos(:,kPast) = measUpt;
            
            for i=so3StateList
                xhatPos.(i)(:,:,k) = xhatPri.(i)(:,:,k)*vec2rot(measUpt(idx.(i)));
                bigphi(idx.(i), idx.(i)) = vec2jac(-measUpt(idx.(i)));
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k) + measUpt(idx.vecState);
            
            if knob.meas.covlim
                H.comb2 = [H.comb; H.covlim];
                R.comb2 = [R.comb' R.covlim];
                K2 = PhatPri(:,:,k)*H.comb2'/(H.comb2*PhatPri(:,:,k)*H.comb2' + diag(R.comb2));
                PhatPos(:,:,k) = bigphi*(I_N-K2*H.comb2)*PhatPri(:,:,k)*bigphi';
            else
                PhatPos(:,:,k) = bigphi*(I_N-K*H.comb)*PhatPri(:,:,k)*bigphi';
            end
        else
            for i=so3StateList
                xhatPos.(i)(:,:,k) = xhatPri.(i)(:,:,k);
            end
            xhatPos.vec(:,k) = xhatPri.vec(:,k);
            PhatPos(:,:,k) = PhatPri(:,:,k);
        end

        %% Constraint update
        N.cstr_k = (knob.cstr.thighlength > 0) * 2 + ...
                   (knob.cstr.hingeknee > 0) * 2;

        n_LT = xhatPos.vec(idx.vec.posPV,k) + xhatPos.W_R_PV(:,:,k)*D.PV_p_LH - ...
                    xhatPos.W_R_LS(:,:,k)*D.LS_p_LK - xhatPos.vec(idx.vec.posLS,k);
        n_RT = xhatPos.vec(idx.vec.posPV,k) + xhatPos.W_R_PV(:,:,k)*D.PV_p_RH - ...
                    xhatPos.W_R_RS(:,:,k)*D.RS_p_RK - xhatPos.vec(idx.vec.posRS,k);

        if knob.cstr.kneerom
            alphaLK = atan2(-dot(n_LT, xhatPos.W_R_LS(:,3,k)), ...
                            -dot(n_LT, xhatPos.W_R_LS(:,1,k))) + 0.5*pi;
            alphaRK = atan2(-dot(n_RT, xhatPos.W_R_RS(:,3,k)), ...
                            -dot(n_RT, xhatPos.W_R_RS(:,1,k))) + 0.5*pi;
            N.cstr_k = N.cstr_k ...
               + ((alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)) ...
               + ((alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax));
        end
 
        D.comb = zeros(N.cstr_k, N.state);
        dy = zeros(N.cstr_k, 1);
        N.cstr_k = 0;
            
        % thigh length constraint
        if knob.cstr.thighlength
            D.comb(N.cstr_k+1,idx.W_R_PV) = -2*n_LT'*xhatPos.W_R_PV(:,:,k)*D.PV_phat_LH;
            D.comb(N.cstr_k+1,idx.W_R_LS) = 2*n_LT'*xhatPos.W_R_LS(:,:,k)*D.LS_phat_LK;
            D.comb(N.cstr_k+1,idx.posPV) = 2*n_LT';
            D.comb(N.cstr_k+1,idx.posLS) = -2*n_LT';
            D.comb(N.cstr_k+2,idx.W_R_PV) = -2*n_RT'*xhatPos.W_R_PV(:,:,k)*D.PV_phat_RH;
            D.comb(N.cstr_k+2,idx.W_R_RS) = 2*n_RT'*xhatPos.W_R_RS(:,:,k)*D.RS_phat_RK;
            D.comb(N.cstr_k+2,idx.posPV) = 2*n_RT';
            D.comb(N.cstr_k+2,idx.posRS) = -2*n_RT';
            dy(N.cstr_k+(1:2),:) = [d.ltl - n_LT'*n_LT; d.rtl - n_RT'*n_RT];
            N.cstr_k = N.cstr_k + 2;
        end
            
        % hinge knee joint constraint
        if knob.cstr.hingeknee
            W_r_LS_y = xhatPos.W_R_LS(:,:,k) * D.i_y;
            W_r_RS_y = xhatPos.W_R_RS(:,:,k) * D.i_y;
            
            D.comb(N.cstr_k+1,idx.W_R_PV) = -W_r_LS_y' * xhatPos.W_R_PV(:,:,k) * D.PV_phat_LH;
            D.comb(N.cstr_k+1,idx.W_R_LS) = W_r_LS_y' * xhatPos.W_R_LS(:,:,k) * D.LS_phat_LK ...
                              - n_LT' * xhatPos.W_R_LS(:,:,k) * D.ihat_y;
            D.comb(N.cstr_k+1,idx.posPV) = W_r_LS_y';
            D.comb(N.cstr_k+1,idx.posLS) = -W_r_LS_y';
            D.comb(N.cstr_k+2,idx.W_R_PV) = -W_r_RS_y' * xhatPos.W_R_PV(:,:,k) * D.PV_phat_RH;
            D.comb(N.cstr_k+2,idx.W_R_RS) = W_r_RS_y' * xhatPos.W_R_RS(:,:,k) * D.RS_phat_RK ...
                              - n_RT' * xhatPos.W_R_RS(:,:,k) * D.ihat_y;
            D.comb(N.cstr_k+2,idx.posPV) = W_r_RS_y';
            D.comb(N.cstr_k+2,idx.posRS) = -W_r_RS_y';

            dy(N.cstr_k+(1:2),:) = [d.lkh - W_r_LS_y'*n_LT; ...
                                    d.rkh - W_r_RS_y'*n_RT];
            N.cstr_k = N.cstr_k + 2;
        end
            
        % knee range of motion
        if knob.cstr.kneerom
            if (alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)
                alphaLK2 = min(max(alphaLK, fOpt.alphaLKmin), fOpt.alphaLKmax);
                a = [-sin(alphaLK2-pi/2); 0; cos(alphaLK2-pi/2)];
                i = xhatPos.W_R_LS(:,:,k) * a;

                D.comb(N.cstr_k+1,idx.W_R_PV) = -i' * xhatPos.W_R_PV(:,:,k) * D.PV_phat_LH;
                D.comb(N.cstr_k+1,idx.W_R_LS) = -n_LT' * xhatPos.W_R_LS(:,:,k) * hat(a) ...
                    + i' * xhatPos.W_R_LS(:,:,k) * D.LS_phat_LK;
                D.comb(N.cstr_k+1,idx.posPV) = i';
                D.comb(N.cstr_k+1,idx.posLS) = -i';
                dy(N.cstr_k+1,:) = d.lkrom - n_LT' * i;

                N.cstr_k = N.cstr_k + 1;
            end
            if (alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax)
                alphaRK2 = min(max(alphaRK, fOpt.alphaRKmin), fOpt.alphaRKmax);
                a = [-sin(alphaRK2-pi/2); 0; cos(alphaRK2-pi/2)];
                i = xhatPos.W_R_RS(:,:,k) * a;

                D.comb(N.cstr_k+1,idx.W_R_PV) = -i' * xhatPos.W_R_PV(:,:,k) * D.PV_phat_RH;
                D.comb(N.cstr_k+1,idx.W_R_RS) = -n_RT' * xhatPos.W_R_RS(:,:,k) * hat(a) ...
                    + i' * xhatPos.W_R_RS(:,:,k) * D.RS_phat_RK;
                D.comb(N.cstr_k+1,idx.posPV) = i';
                D.comb(N.cstr_k+1,idx.posRS) = -i';
                dy(N.cstr_k+1,:) = d.rkrom - n_RT' * i;

                N.cstr_k = N.cstr_k + 1;
            end
        end
        
        % apply constraint update
        if N.cstr_k > 0    
            K = PhatPos(:,:,k)*D.comb'/(D.comb*PhatPos(:,:,k)*D.comb');
            if knob.cstr.Pupdate
                Ptilde(:,:,k) = (I_N-K*D.comb)*PhatPos(:,:,k);
            else
                Ptilde(:,:,k) = PhatPos(:,:,k);
            end

            measUpt = K*dy;
            debug_dat.measUptTilde(:,kPast) = measUpt;
                        
            for i=so3StateList
                xtilde.(i)(:,:,k) = xhatPos.(i)(:,:,k)*vec2rot(measUpt(idx.(i)));
            end
            xtilde.vec(:,k) = xhatPos.vec(:,k) + measUpt(idx.vecState);
        else
            for i=so3StateList
                xtilde.(i)(:,:,k) = xhatPos.(i)(:,:,k);
            end
            xtilde.vec(:,k) = xhatPos.vec(:,k);
            Ptilde(:,:,k) = PhatPos(:,:,k);
        end
    end
    
    %% remove offset state (k=0)
    for i=so3StateList
        xhatPri.(i) = xhatPri.(i)(:,:,2:end);
        xhatPos.(i) = xhatPos.(i)(:,:,2:end);
        xtilde.(i) = xtilde.(i)(:,:,2:end);
    end
    xhatPri.vec = xhatPri.vec(:,2:end);
    xhatPos.vec = xhatPos.vec(:,2:end);
    xtilde.vec = xtilde.vec(:,2:end);
    PhatPri = PhatPri(:,:,2:end);
    PhatPos = PhatPos(:,:,2:end);
    Ptilde = Ptilde(:,:,2:end);
    
    %% compute debug information
    LTIBz = squeeze(xtilde.W_R_LS(:,3,:))';
    RTIBz = squeeze(xtilde.W_R_RS(:,3,:))';
    PELVy = squeeze(xtilde.W_R_PV(:,2,:))';
    debug_dat.LFEO = xtilde.vec(idx.vec.posLS,:)' + body.LS_d * LTIBz;
    debug_dat.RFEO = xtilde.vec(idx.vec.posRS,:)' + body.RS_d * RTIBz;
    debug_dat.LFEP = xtilde.vec(idx.vec.posPV,:)' + body.PV_d/2 * PELVy;
    debug_dat.RFEP = xtilde.vec(idx.vec.posPV,:)' - body.PV_d/2 * PELVy;
    
    R_LFEM = zeros(3,3,N.samples);
    R_RFEM = zeros(3,3,N.samples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(xtilde.W_R_LS(:,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(xtilde.W_R_RS(:,2,:));
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