%> EKF using Lie group/algebra representation
%> 3 KMUs presumably worn on the body in the following configuration: 
%> mid pelvis, left ankle, right ankle
%>
%> Author: Luke Wicent Sy
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

function [ xhat_pri, xtilde, debug_dat ] = lieekf_3_kmus_v3(x0, P0, ...
    bodyAccMP, bIsStatMP, qMP, wbodyMP, ...
    bodyAccLA, bIsStatLA, qLA, wbodyLA, ...
    bodyAccRA, bIsStatRA, qRA, wbodyRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)
    [nSamples, ~] = size(bodyAccMP);
    
    fOpt = struct('fs', 60);
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    % step = struct('RPV', bIsStatMP, 'LSK', bIsStatLA, 'RSK', bIsStatRA);
    qOri = struct('W_RPV_T', qMP, 'W_LSK_T', qLA, 'W_RSK_T', qRA);
    % wOri = struct('RPV', wMP, 'LSK', wLA, 'RSK', wRA);
    u = [bodyAccMP'; bodyAccLA'; bodyAccRA';
         wbodyMP'; wbodyLA'; wbodyRA'];
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
%     xhat_pos.W_RPV_T = zeros(4,4,nSamples);
%     xhat_pos.W_LSK_T = zeros(4,4,nSamples);
%     xhat_pos.W_RSK_T = zeros(4,4,nSamples);
%     xhat_pos.W_RPV_T(1:3,1:3,:) = quat2rotm(qMP); 
%     xhat_pos.W_RPV_T(4,4,:) = 1;
%     xhat_pos.W_LSK_T(1:3,1:3,:) = quat2rotm(qLA); 
%     xhat_pos.W_LSK_T(4,4,:) = 1;
%     xhat_pos.W_RSK_T(1:3,1:3,:) = quat2rotm(qRA); 
%     xhat_pos.W_RSK_T(4,4,:) = 1;

    %% State initialization
    xhat_pri = struct();
    xhat_pos = struct();
    xtilde = struct();
    
    state1List = {'W_RPV_T', 'W_LSK_T', 'W_RSK_T'};
    state1ListN = length(state1List);
    
    for i=1:3
        sname = state1List{i};
        xhat_pri.(sname) = nan(4,4,nSamples+1);
        xhat_pos.(sname) = nan(4,4,nSamples+1);
        xtilde.(sname) = nan(4,4,nSamples+1);
    end
    xhat_pri.v = nan(18,nSamples+1);
    xhat_pos.v = nan(18,nSamples+1);
    xtilde.v = nan(18,nSamples+1);
    
    Fv = [eye(9,9) zeros(9,9);
         zeros(9,18)];
    Gv = [dt*eye(9,9) zeros(9,9);
         zeros(9,9) eye(9,9)];
     
    %% Set n=1 states (initial state taken from input)
    % se(3) pos and ori states
    x0Offset = [0, 10, 20];
    for i=1:state1ListN
        sname = state1List{i};
        offset = x0Offset(i);
        
        xtilde.(sname)(:,:,1) = zeros(4,4);
        xtilde.(sname)(1:3,1:3,1) = quat2rotm(x0(offset+(7:10))');
        xtilde.(sname)(1:3,4,1) = x0(offset+(1:3));
        xtilde.(sname)(4,4,1) = 1;
    end
    % vel and ang. vel
    xtilde.v(:,1) = [x0([4:6 14:16 24:26]); zeros(9,1)];
    vList = {[1:3 10:12] [4:6 13:15] [7:9 16:18]};
    
    % Iteration
    xiMatDebug = zeros(4,4,nSamples);
    for n=2:(nSamples+1)
        %% Prediction
        for i=1:state1ListN
            sname = state1List{i};
            xi = xtilde.v(vList{i},n-1);
            xi(1:3) = quatrotate(qOri.(sname)(n-1,:), xi(1:3)')';
            xiMatDebug(:,:,n-1) = expm(liese3.caret(xi*dt));
            xhat_pri.(sname)(:,:,n) = xtilde.(sname)(:,:,n-1)*xiMatDebug(:,:,n-1);
%             xhat_pri.(sname)(1:3,1:3,n) = quat2rotm(qOri.(sname)(n-1,:));
        end
        xhat_pri.v(:,n) = Fv*xtilde.v(:,n-1) + Gv*u(:,n-1);
        
%         pRPV = pRPV + vRPV*dt + gfrAccMP(n,:)*dt2;
%         pLSK = pLSK + vLSK*dt + gfrAccLA(n,:)*dt2;
%         pRSK = pRSK + vRSK*dt + gfrAccRA(n,:)*dt2;
%         
%         vRPV = vRPV + gfrAccMP(n,:)*dt;
%         vLSK = vLSK + gfrAccLA(n,:)*dt;
%         vRSK = vRSK + gfrAccRA(n,:)*dt;
%         
%         xhat_pos.W_RPV_T(1:3,4,n) = pRPV';
%         xhat_pos.W_LSK_T(1:3,4,n) = pLSK';
%         xhat_pos.W_RSK_T(1:3,4,n) = pRSK';
        for i=1:state1ListN
            sname = state1List{i};
            xtilde.(sname)(:,:,n) = xhat_pri.(sname)(:,:,n);
        end
        xtilde.v(:,n) = xhat_pri.v(:,n);
    end
    
    for i=1:state1ListN
        sname = state1List{i};
        xhat_pri.(sname) = xhat_pri.(sname)(:,:,2:end);
        xhat_pos.(sname) = xhat_pos.(sname)(:,:,2:end);
        xtilde.(sname) = xtilde.(sname)(:,:,2:end);
    end
    xhat_pri.v = xhat_pri.v(:,2:end);
    xhat_pos.v = xhat_pos.v(:,2:end);
    xtilde.v = xtilde.v(:,2:end);
    
    LTIBz = squeeze(xtilde.W_LSK_T(1:3,3,:))';
    RTIBz = squeeze(xtilde.W_RSK_T(1:3,3,:))';
    PELVy = squeeze(xtilde.W_RPV_T(1:3,2,:))';
    debug_dat.LFEO = squeeze(xtilde.W_LSK_T(1:3,4,:))' + dLTibia * LTIBz;
    debug_dat.RFEO = squeeze(xtilde.W_RSK_T(1:3,4,:))' + dRTibia * RTIBz;
    debug_dat.LFEP = squeeze(xtilde.W_RPV_T(1:3,4,:))' + dPelvis/2 * PELVy;
    debug_dat.RFEP = squeeze(xtilde.W_RPV_T(1:3,4,:))' - dPelvis/2 * PELVy;
    
    R_LFEM = zeros(3,3,nSamples);
    R_RFEM = zeros(3,3,nSamples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(xtilde.W_LSK_T(1:3,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(xtilde.W_RSK_T(1:3,2,:));
    RFEM_x = cross(RFEM_y, RFEM_z);
    R_RFEM(:,3,:) = RFEM_z ./ vecnorm(RFEM_z, 2, 1);
    R_RFEM(:,2,:) = RFEM_y ./ vecnorm(RFEM_y, 2, 1);
    R_RFEM(:,1,:) = RFEM_x ./ vecnorm(RFEM_x, 2, 1);
    debug_dat.qLTH = rotm2quat(R_LFEM);
    debug_dat.qRTH = rotm2quat(R_RFEM);
end