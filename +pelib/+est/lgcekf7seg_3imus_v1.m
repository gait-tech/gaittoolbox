function [ xtilde, debug_dat ] = lgcekf7seg_3imus_v1(x0, P0, ...
    B_a_, step, W_R_, B_w_, body, uwb_mea, options)
% Constrained EKF using Lie group/algebra representation
% Representation: SE(3)^3 + R^9
% 3 IMUs presumably worn on the body in the following configuration: 
% mid pelvis, left ankle, right ankle
% 
% Coding notation based on <a href="http://paulfurgale.info/news/2014/6/9/representing-robot-pose-the-good-the-bad-and-the-ugly">link</a>
%
% More detail about options
%      fs: sampling frequency of the magnetic and inertial measurement units
%      applyPred: 3 digit ZYX
%      applyMeas: 3 digit ZYX
%      applyCstr: 3 digit ZYX
%
% :param x0: initial state in the GFR
% :param P0: initial covariance
% :param B_a_: acceleration of PV, LS, RS in the body frame
% :param step: boolean vector of PV, LS, RS indicating step detection
% :param W_R_: PV, LS, RS orientation in the GFR (rotm)
% :param B_w_: PV, LS, RS angular velocity in the body frame
% :param body: Length of PV_d (pelvis), RT_d and LT_d (r/l femur), RS_d and LS_d (r/l tibia)
% :param uwb_mea: a structure containing the range measurements (m) between
% :param options: struct containing the estimator settings:
%
% .. Author: - Luke Wicent Sy (GSBME, 2019 Oct 31)
    N = {};
    [N.samples, ~] = size(B_a_.PV);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'applyPred', 1, 'applyMeas', 1, 'applyCstr', 1 );
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% Feature on off
    knob = struct();
       
    %% State and error covariance initialization
    se3StateList = ["W_T_PV", "W_T_LF", "W_T_RF"];
    bodyList = ["PV", "LF", "RF"];
    N.se3State = length(se3StateList)*3;    
    N.r3State = 3*3;
    N.state = N.se3State + N.r3State;
    
    % calculate indices
    idx = struct('se3State', 1:9, ...
                 'W_T_PV', 1:3, 'W_T_LF', 4:6, 'W_T_RF', 7:9, ...
                 'vec', struct() );
%     for i=["pos", "vel", "avel"]
%         for j=bodyList
%             k = sprintf('%s%s', i, j);
%             idx.(k) = (1:3) + N.state;
%             idx.vec.(k) = (1:3) + N.state - N.se3State;
%             N.state = N.state + 3;
%         end
%     end
%     idx.vecState = (N.se3State+1):N.state;
    
    g = [0 0 9.81]';        % gravity
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    I_N = eye(N.state);
    
    xhatPri = struct();     % prediction update state
    xhatPos = struct();     % measurement update state
    xtilde = struct();      % constraint update state
    for i=se3StateList
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
%     xtilde.W_T_PV(:,:,1) = quat2rotm(x0(07:10)');
%     xtilde.W_T_LF(:,:,1) = quat2rotm(x0(17:20)');
%     xtilde.W_T_RF(:,:,1) = quat2rotm(x0(27:30)');
%     xtilde.vec(:,1) = [x0([01:03 11:13 21:23]); ... % position
%                        x0([04:06 14:16 24:26]); ... % velocity
%                        zeros(9,1)]; % angular velocity
    if isscalar(P0)
        P0 = I_N*P0;
    end
    Ptilde(:,:,1) = P0;
    
    %% Prediction step initialization
    
    %% Measurement step initialization
    
    %% Constraint step initialization
    D = struct('i_x', [1 0 0]', 'i_y', [0 1 0]', 'i_z', [0 0 1]');
    D.ihat_x = hat(D.i_x); D.ihat_y = hat(D.i_y); D.ihat_z = hat(D.i_z); 
    d = {};
        
    %% Iteration
    for k=2:(N.samples+1)
        kPast = k-1;
        
        %% Prediction update
%         if knob.pred.useStateinWFrameConv
%             u(:,kPast) = [xtilde.W_T_PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
%                  xtilde.W_T_LF(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
%                  xtilde.W_T_RF(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
%                  xtilde.W_T_PV(:,:,kPast)*B_w2_.PV(kPast,:)'; ...
%                  xtilde.W_T_LF(:,:,kPast)*B_w2_.LS(kPast,:)'; ...
%                  xtilde.W_T_RF(:,:,kPast)*B_w2_.RS(kPast,:)'];
%         else
%             u(:,kPast) = [W_R_.PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
%                  W_R_.LS(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
%                  W_R_.RS(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
%                  W_R_.PV(:,:,kPast)*B_w2_.PV(kPast,:)'; ...
%                  W_R_.LS(:,:,kPast)*B_w2_.LS(kPast,:)'; ...
%                  W_R_.RS(:,:,kPast)*B_w2_.RS(kPast,:)'];
%         end
%         
%         F.AdG = eye(N.state, N.state);
%         bigphi = eye(N.state, N.state);
%         for i=se3StateList
%             j = string(i{1}(end-1:end));
%             
%             if knob.pred.ori
%                 xi = B_w2_.(j)(kPast,:)'; % B_avel
%             else
%                 xi = zeros(3,1);
%             end
%             
%             xhatPri.(i)(:,:,k) = xtilde.(i)(:,:,kPast)*vec2rot(xi*dt);
%             
%             bigphi(idx.(i), idx.(i)) = vec2jac(-xi*dt);
%             % for SO(3) Ad(X) = X;
%             % F.AdG(idx.(i), idx.(i)) = tranAd(vec2rot(-xi*dt));
%             F.AdG(idx.(i), idx.(i)) = vec2rot(-xi*dt);
%         end
%         
%         if knob.pred.posvelangvel
%             xhatPri.vec(:,k) = F.vec*xtilde.vec(:,kPast) + G.vec*u(:,kPast);
%         else
%             xhatPri.vec(:,k) = xtilde.vec(:,kPast);
%         end
%         
%         F.comb = F.AdG + bigphi*F.curlyC;
%         PhatPri(:,:,k) = F.comb*Ptilde(:,:,kPast)*F.comb' + bigphi*Q.comb*bigphi';
        
        %% Measurement update

        %% Constraint update

    end
    
    %% remove offset state (k=0)
    for i=se3StateList
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
%     debug_dat.u = u;
    debug_dat.xhatPri = xhatPri;
    debug_dat.xhatPos = xhatPos;
    debug_dat.xtilde = xtilde;
    debug_dat.PhatPri = PhatPri;
    debug_dat.PhatPos = PhatPos;
    debug_dat.Ptilde = Ptilde;
end