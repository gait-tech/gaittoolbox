%% Reconstruct motion using kf_3_kmus_v2 (KF with constraints)
% Experiment 1b:
% Hinge Knee Joint Constraint @ Vicon Frame.
% Input vicon generated IMU data + gradual increase of noise
% Sample Run:
% dataV = loadBVHasAnimatedata('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'mm');
% dataS = loadSensors('Data/s1/Acting1_Xsens_AuxFields.sensors');
% animateTrialSample05(dataV, dataS, 60);
% animateTrialSample05('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'Data/s1/Acting1_Xsens_AuxFields.sensors', 60);


% function results = tcdExperiment01b(fnameV, fnameS)
    %% Inputs and Input Check
    fnameV = 'totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh';
    fnameS = 'totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors';
    fnameCIB = 'totalcapture/imu/s1/s1_acting1_calib_imu_bone.txt';
    fnameCIR = 'totalcapture/imu/s1/s1_acting1_calib_imu_ref.txt';
    % Check function input
    validateattributes(fnameV, {'char', 'struct'}, {});
    validateattributes(fnameS, {'char', 'struct'}, {});
    
    %% Initialization
    % Initialize variables and libraries
    fs = 60;
    
    % Load calibration data
    calibIB = mocapdb.XsensBody.loadCalib(fnameCIB);
    calibIR = mocapdb.XsensBody.loadCalib(fnameCIR);
    
    % Load video orientation and position for each body segment
    if ischar(fnameV)
        dataV = mocapdb.BVHBody.loadBVHFile(fnameV, 'mm');
    else
        dataV = fnameV;
    end     
    
    % Load sensor orientation and position for each body segment
    if ischar(fnameS)
        dataS = mocapdb.XsensBody.loadSensorFile(fnameS);
    else
        dataS = fnameS;
    end
    
    n = min(dataV.nSamples, dataS.nSamples);
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
        'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
    m = containers.Map(key, val);
    for i=1:length(key)
        dataS.(key{i}) = dataS.(key{i})(1:n,:);
        dataV.(val{i}) = dataV.(val{i})(1:n,:)/1000;
%         dataV.(val{i}) = quatrotate(calibIR.(key{i}).ori, ...
%             dataV.(val{i})(1:n,:)/1000);
%         v2 = strcat('q', val{i});
%         dataV.(v2) = quatmultiply(quatconj(calibIR.(key{i}).ori), ...
%             dataV.(v2)(1:n,:)/1000);
    end
    dataV.posUnit = 'm';
    
    MIDPEL_act = [mean([dataV.LeftUpLeg(:,1) dataV.RightUpLeg(:,1)], 2),...
                  mean([dataV.LeftUpLeg(:,2) dataV.RightUpLeg(:,2)], 2),...
                  mean([dataV.LeftUpLeg(:,3) dataV.RightUpLeg(:,3)], 2)];
%     gfr_vel_MP_act = [gradient(MIDPEL_act(:,1), 1/fs), ...
%                       gradient(MIDPEL_act(:,2), 1/fs), ...
%                       gradient(MIDPEL_act(:,3), 1/fs)];
%     gfr_acc_MP_act = [gradient(gfr_vel_MP_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_MP_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_MP_act(:,3), 1/fs)];
%     gfr_vel_LA_act = [gradient(dataV.LeftFoot(:,1), 1/fs), ...
%                       gradient(dataV.LeftFoot(:,2), 1/fs), ...
%                       gradient(dataV.LeftFoot(:,3), 1/fs)];
%     gfr_acc_LA_act = [gradient(gfr_vel_LA_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_LA_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_LA_act(:,3), 1/fs)];
%     gfr_vel_RA_act = [gradient(dataV.RightFoot(:,1), 1/fs), ...
%                       gradient(dataV.RightFoot(:,2), 1/fs), ...
%                       gradient(dataV.RightFoot(:,3), 1/fs)];
%     gfr_acc_RA_act = [gradient(gfr_vel_RA_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_RA_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_RA_act(:,3), 1/fs)];
    gfr_vel_MP_act = diff(MIDPEL_act, 1, 1)*fs;
    gfr_vel_MP_act = [gfr_vel_MP_act(1,:); gfr_vel_MP_act];
    gfr_acc_MP_act = [0 0 0; diff(MIDPEL_act, 2, 1)*fs*fs];
    gfr_vel_LA_act = diff(dataV.LeftFoot, 1, 1)*fs;
    gfr_vel_LA_act = [gfr_vel_LA_act(1,:); gfr_vel_LA_act];
    gfr_acc_LA_act = [0 0 0; diff(dataV.LeftFoot, 2, 1)*fs*fs];
    gfr_vel_RA_act = diff(dataV.RightFoot, 1, 1)*fs;
    gfr_vel_RA_act = [gfr_vel_RA_act(1,:); gfr_vel_RA_act];
    gfr_acc_RA_act = [0 0 0; diff(dataV.RightFoot, 2, 1)*fs*fs];
    
    x0_pos_MP = MIDPEL_act(1,:);
    x0_pos_LA = dataV.LeftFoot(1,:);
    x0_pos_RA = dataV.RightFoot(1,:);
    x0_vel_MP = gfr_vel_MP_act(1,:);
    x0_vel_LA = gfr_vel_LA_act(1,:);
    x0_vel_RA = gfr_vel_RA_act(1,:);
    sigma_acc = 0.5;
    
    WIN_SECS = 0.25;
    VAR_WIN  = floor(fs*WIN_SECS); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;

    movVarAcc_pelvis = movingvar(sqrt( sum(dataS.Pelvis.acc.^2,2)) ,VAR_WIN);
    bIsStatMP = movVarAcc_pelvis < 0;
    movVarAcc_lankle = movingvar(sqrt( sum(dataS.L_LowLeg.acc.^2,2)) ,VAR_WIN);
    bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(dataS.R_LowLeg.acc.^2,2)) ,VAR_WIN);
    bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
    
    movVarAcc_pelvis_act = movingvar(sqrt( sum(gfr_acc_MP_act.^2,2)) ,VAR_WIN);
    bIsStatMP_act = movVarAcc_pelvis_act < 0;
    movVarAcc_lankle_act = movingvar(sqrt( sum(gfr_acc_LA_act.^2,2)) ,VAR_WIN);
    bIsStatLA_act = movVarAcc_lankle_act < ACC_VAR_THRESH;
    movVarAcc_rankle_act = movingvar(sqrt( sum(gfr_acc_RA_act.^2,2)) ,VAR_WIN);
    bIsStatRA_act = movVarAcc_rankle_act < ACC_VAR_THRESH;
    
    % orientation of body in Vicon frame
    qPelvisEst0 = quatmultiply((calibIR.Pelvis.ori), ...
        quatmultiply(dataS.Pelvis.ori, quatconj(calibIB.Pelvis.ori)));
    qLankleEst0 = quatmultiply((calibIR.L_LowLeg.ori), ...
        quatmultiply(dataS.L_LowLeg.ori, quatconj(calibIB.L_LowLeg.ori)));
    qRankleEst0 = quatmultiply((calibIR.R_LowLeg.ori), ...
        quatmultiply(dataS.R_LowLeg.ori, quatconj(calibIB.R_LowLeg.ori)));
    % orientation of body in World frame
%     qPelvisEst = quatmultiply(dataS.Pelvis.ori, ...
%                               quatconj(calibIB.Pelvis.ori));
%     qLankleEst = quatmultiply(dataS.L_LowLeg.ori, ...
%                               quatconj(calibIB.L_LowLeg.ori));
%     qRankleEst = quatmultiply(dataS.R_LowLeg.ori, ...
%                               quatconj(calibIB.R_LowLeg.ori));

    % Align body frame convention to biomechanics convention
    % FROM: x - right, y - back, z - down
    % TO: x - front, y - left, z - up
    qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
    qPelvisEst = quatmultiply(qPelvisEst0, qTCD2BM);
    qLankleEst = quatmultiply(qLankleEst0, qTCD2BM);
    qRankleEst = quatmultiply(qRankleEst0, qTCD2BM);
    qPelvisAct = quatmultiply(dataV.qHips, qTCD2BM);
    qLankleAct = quatmultiply(dataV.qLeftLeg, qTCD2BM);
    qRankleAct = quatmultiply(dataV.qRightLeg, qTCD2BM);
    
    gfr_acc_MP = quatrotate(quatconj(dataS.Pelvis.ori), ...
                            dataS.Pelvis.acc) - [0 0 9.81];
    gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
    gfr_acc_LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
                            dataS.L_LowLeg.acc) - [0 0 9.81];
    gfr_acc_LA = quatrotate(quatconj(calibIR.L_LowLeg.ori), gfr_acc_LA);
    gfr_acc_RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
                            dataS.R_LowLeg.acc) - [0 0 9.81];
    gfr_acc_RA = quatrotate(quatconj(calibIR.R_LowLeg.ori), gfr_acc_RA);
    
    fc = 10;
    [lpf_b, lpf_a] = butter(6, fc/(fs/2));
    gfr_acc_MP_filt = filter(lpf_b, lpf_a, gfr_acc_MP);
    gfr_acc_LA_filt = filter(lpf_b, lpf_a, gfr_acc_LA);
    gfr_acc_RA_filt = filter(lpf_b, lpf_a, gfr_acc_RA);
    
    %% -----------------------------------------------------------------------
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_act-dataV.LeftFoot), 2, 2);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_act-dataV.RightFoot), 2, 2);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.RightFoot-dataV.LeftFoot), 2, 2);

    d_pelvis = norm(dataV.RightUpLeg(1,:) - dataV.LeftUpLeg(1,:));
    d_rfemur = norm(dataV.RightUpLeg(1,:) - dataV.RightLeg(1,:));
    d_lfemur = norm(dataV.LeftUpLeg(1,:) - dataV.LeftLeg(1,:));
    d_rtibia = norm(dataV.RightLeg(1,:) - dataV.RightFoot(1,:));
    d_ltibia = norm(dataV.LeftLeg(1,:) - dataV.LeftFoot(1,:));
    
    P = eye(18, 18)*100;
    Q_acc = sigma_acc*100;
    
    % validation initialization
    idx = 1:240;
    actBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', ...
                         'xyzColor', {'m', 'y', 'c'}}); 

    N_TRIALS = 1;
    resultsIdx = 1;
    for i=1.5:1.5
        for j=1:N_TRIALS
            noiseLen = length(gfr_acc_MP_act(:,1));
            gfr_acc_MP_actnoise = gfr_acc_MP_act + ...
                i*[randn(noiseLen,1) randn(noiseLen,1) randn(noiseLen,1)];
            gfr_acc_LA_actnoise = gfr_acc_LA_act + ...
                i*[randn(noiseLen,1) randn(noiseLen,1) randn(noiseLen,1)];
            gfr_acc_RA_actnoise = gfr_acc_RA_act + ...
                i*[randn(noiseLen,1) randn(noiseLen,1) randn(noiseLen,1)];

            % No constraint, Vicon acc, Vicon ori
            [ x_pri_v2, x_pos_v2, t_dat_v2 ] = kf_3_kmus_v2(fs, ...
                Q_acc, Q_acc, Q_acc, P, ...
                x0_pos_MP, x0_vel_MP, gfr_acc_MP_actnoise, ...
                bIsStatMP_act, qPelvisAct(2:end,:), ...
                x0_pos_LA, x0_vel_LA, gfr_acc_LA_actnoise, ...
                bIsStatLA_act, qLankleAct(2:end,:), ...
                x0_pos_RA, x0_vel_RA, gfr_acc_RA_actnoise, ...
                bIsStatRA_act, qRankleAct(2:end,:), ...
                d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
                true, false, true, false, false);

            estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                       'lnSymbol', '--', 'ptSymbol', 'o', ...
                       'xyzColor', {'r', 'g', 'b'}, ...
                       'MIDPEL', x_pos_v2(idx,1:3), ...
                       'LFEP', t_dat_v2.LFEP(idx,:), ...
                       'LFEO', t_dat_v2.LFEO(idx,:), ...
                       'LTIO', x_pos_v2(idx,7:9), ...
                       'RFEP', t_dat_v2.RFEP(idx,:), ...
                       'RFEO', t_dat_v2.RFEO(idx,:), ...
                       'RTIO', x_pos_v2(idx,13:15), ...
                       'qRPV', qPelvisAct(idx+1,:), ...
                       'qRTH', t_dat_v2.qRTH(idx,:), ...
                       'qLTH', t_dat_v2.qLTH(idx,:), ...
                       'qRSK', qRankleAct(idx+1,:), ...
                       'qLSK', qLankleAct(idx+1,:));
            trial(j) = estBody.diffRMSE(actBody);
        end
        
%         struct2csvstr(trial)
        
        results(resultsIdx) = pelib.diffRMSEMean(trial);
        resultsIdx = resultsIdx + 1;
    end
    
    struct2csvstr(results, true)

    actState = [MIDPEL_act gfr_vel_MP_act ...
                dataV.LeftFoot gfr_vel_LA_act ...
                dataV.RightFoot gfr_vel_RA_act];

    %% --------------------------------------------------------------------
    %  Further Validation
    % Static Plots
    updateFigureContents('Position');
    pelib.viz.plotPosition({estBody, actBody}, {'MIDPEL', 'LTIO', 'RTIO'});
    
    updateFigureContents('Animation Freeze');
    grid; view(0, 90); hold on;
    for i=idx(1):20:idx(end)
        pelib.viz.plotLowerBody(estBody, i);
        pelib.viz.plotLowerBody(actBody, i);
    end
    
    % Animation
%     updateFigureContents('Animation');
%     xlabel('x'); ylabel('y'); zlabel('z');
%     estBodyLimits = [estBody.xlim() estBody.ylim() estBody.zlim()];
%     for i=idx
%         clf; grid;
%         xlim(estBodyLimits(1:2)); 
%         ylim(estBodyLimits(3:4)); 
%         zlim(estBodyLimits(5:6));  
%         view(0, 180);
%         pelib.viz.plotLowerBody(estBody, i);
%         pause(1/1000);
%     end

%     updateFigureContents('Animation');
%     xlabel('x'); ylabel('y'); zlabel('z');
%     actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()];
%     for i=idx
%         clf; grid; 
%         xlim(actBodyLimits(1:2)); 
%         ylim(actBodyLimits(3:4)); 
%         zlim(actBodyLimits(5:6));  
%         view(0, 180);
%         pelib.viz.plotLowerBody(actBody, i);
%         pause(1/1000);
%     end
% end