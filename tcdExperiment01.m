%% Reconstruct motion using kf_3_kmus (KF with constraints)
% Steps:
% - calculate orientation (from Xsens or CAHRS)
% - calculate root position by double integration
% - measure angle and position RMSE
% - reconstruct skeleton to feed into visualization
% Sample Run:
% dataV = loadBVHasAnimatedata('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'mm');
% dataS = loadSensors('Data/s1/Acting1_Xsens_AuxFields.sensors');
% animateTrialSample05(dataV, dataS, 60);
% animateTrialSample05('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'Data/s1/Acting1_Xsens_AuxFields.sensors', 60);

% function results = tcdExperiment01(fnameV, fnameS)
    %% Inputs and Input Check
    fnameV = 'TotalCaptureDataset/s1/acting1_BlenderZXY_YmZ.bvh';
    fnameS = 'TotalCaptureDataset/s1/Acting1_Xsens_AuxFields.sensors';
    fnameCIB = 'TotalCaptureDataset/IMUCalibration/calib_imu_bone.txt';
    fnameCIR = 'TotalCaptureDataset/IMUCalibration/calib_imu_ref.txt';
%     fnameCIB = 'totalcapture/imu/s1/s1_acting1_calib_imu_bone.txt';
%     fnameCIR = 'totalcapture/imu/s1/s1_acting1_calib_imu_ref.txt';
    % Check function input
    validateattributes(fnameV, {'char', 'struct'}, {});
    validateattributes(fnameS, {'char', 'struct'}, {});
    
    %% Initialization
    % Initialize variables and libraries
    fs = 60;
    
    % Load calibration data
    calibIB = tcdlib.XsensBody.loadCalib(fnameCIB);
    calibIR = tcdlib.XsensBody.loadCalib(fnameCIR);
    
    % Load video orientation and position for each body segment
    if ischar(fnameV)
        dataV = tcdlib.BVHBody.loadBVHFile(fnameV, 'mm');
    else
        dataV = fnameV;
    end     
    
    % Load sensor orientation and position for each body segment
    if ischar(fnameS)
        dataS = tcdlib.XsensBody.loadSensorFile(fnameS);
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
    
    x0_pos_MP = dataV.Hips(1,:);
    x0_pos_LA = dataV.LeftFoot(1,:);
    x0_pos_RA = dataV.RightFoot(1,:);
    x0_vel_MP = [0 0 0];
    x0_vel_LA = [0 0 0];
    x0_vel_RA = [0 0 0];
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
    
    % orientation of body in Vicon frame
    qPelvisEst = quatmultiply((calibIR.Pelvis.ori), ...
        quatmultiply(dataS.Pelvis.ori, quatconj(calibIB.Pelvis.ori)));
    qLankleEst = quatmultiply((calibIR.L_LowLeg.ori), ...
        quatmultiply(dataS.L_LowLeg.ori, quatconj(calibIB.L_LowLeg.ori)));
    qRankleEst = quatmultiply((calibIR.R_LowLeg.ori), ...
        quatmultiply(dataS.R_LowLeg.ori, quatconj(calibIB.R_LowLeg.ori)));
    % orientation of body in IMU frame
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
    qPelvisEst = quatmultiply(qPelvisEst, qTCD2BM);
    qLankleEst = quatmultiply(qLankleEst, qTCD2BM);
    qRankleEst = quatmultiply(qRankleEst, qTCD2BM);
    
    gfr_acc_MP = quatrotate(quatconj(dataS.Pelvis.ori), ...
                            dataS.Pelvis.acc) - [0 0 9.81];
    gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
    gfr_acc_LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
                            dataS.L_LowLeg.acc) - [0 0 9.81];
    gfr_acc_LA = quatrotate(quatconj(calibIR.L_LowLeg.ori), gfr_acc_LA);
    gfr_acc_RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
                            dataS.R_LowLeg.acc) - [0 0 9.81];
    gfr_acc_RA = quatrotate(quatconj(calibIR.R_LowLeg.ori), gfr_acc_RA);
    
    %% -----------------------------------------------------------------------
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    uwb_mea.left_tibia_mid_pelvis = vecnorm((dataV.Hips-dataV.LeftFoot)', 2, 2);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((dataV.Hips-dataV.RightFoot)', 2, 2);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.RightFoot-dataV.LeftFoot)', 2, 2);

    d_pelvis = norm(dataV.RightUpLeg(1,:) - dataV.LeftUpLeg(1,:));
    d_rfemur = norm(dataV.RightUpLeg(1,:) - dataV.RightLeg(1,:));
    d_lfemur = norm(dataV.LeftUpLeg(1,:) - dataV.LeftLeg(1,:));
    d_rtibia = norm(dataV.RightLeg(1,:) - dataV.RightFoot(1,:));
    d_ltibia = norm(dataV.LeftLeg(1,:) - dataV.LeftFoot(1,:));
    
    [ x_pri_v2, x_pos_v2, t_dat_v2 ] = kf_3_kmus_v2(fs, ...
        sigma_acc, sigma_acc, sigma_acc, false, ...
        x0_pos_MP, x0_vel_MP, gfr_acc_MP, bIsStatMP, qPelvisEst, ...
        x0_pos_LA, x0_vel_LA, gfr_acc_LA, bIsStatLA, qLankleEst, ...
        x0_pos_RA, x0_vel_RA, gfr_acc_RA, bIsStatRA, qRankleEst, ...
        d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
        true, false, true, false, false);
    
    %% --------------------------------------------------------------------
    %  Validation
    idx = 1:240;
    estBody = grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
               'lnSymbol', '--', 'ptSymbol', 'o', ...
               'xyzColor', {'r', 'g', 'b'}, ...
               'MIDPEL', x_pos_v2(idx,1:3), ...
               'LFEP', t_dat_v2.LFEP(idx,:), ...
               'LFEO', t_dat_v2.LFEO(idx,:), ...
               'LTIO', x_pos_v2(idx,7:9), ...
               'RFEP', t_dat_v2.RFEP(idx,:), ...
               'RFEO', t_dat_v2.RFEO(idx,:), ...
               'RTIO', x_pos_v2(idx,13:15), ...
               'qRPV', qPelvisEst(idx,:), ...
               'qRTH', t_dat_v2.qRTH(idx,:), ...
               'qLTH', t_dat_v2.qLTH(idx,:), ...
               'qRSK', qRankleEst(idx,:), ...
               'qLSK', qLankleEst(idx,:));
           
    actBody = dataV.togrBody(idx, {'name', 'act', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '.', ...
                             'xyzColor', {'m', 'y', 'c'}});

    gfr_vel_MP_act = [gradient(dataV.Hips(:,1), 1/fs), ...
                      gradient(dataV.Hips(:,2), 1/fs), ...
                      gradient(dataV.Hips(:,3), 1/fs)];
    gfr_acc_MP_act = [gradient(gfr_vel_MP_act(:,1), 1/fs), ...
                      gradient(gfr_vel_MP_act(:,2), 1/fs), ...
                      gradient(gfr_vel_MP_act(:,3), 1/fs)];
    gfr_vel_LA_act = [gradient(dataV.LeftFoot(:,1), 1/fs), ...
                      gradient(dataV.LeftFoot(:,2), 1/fs), ...
                      gradient(dataV.LeftFoot(:,3), 1/fs)];
    gfr_acc_LA_act = [gradient(gfr_vel_LA_act(:,1), 1/fs), ...
                      gradient(gfr_vel_LA_act(:,2), 1/fs), ...
                      gradient(gfr_vel_LA_act(:,3), 1/fs)];
    gfr_vel_RA_act = [gradient(dataV.RightFoot(:,1), 1/fs), ...
                      gradient(dataV.RightFoot(:,2), 1/fs), ...
                      gradient(dataV.RightFoot(:,3), 1/fs)];
    gfr_acc_RA_act = [gradient(gfr_vel_RA_act(:,1), 1/fs), ...
                      gradient(gfr_vel_RA_act(:,2), 1/fs), ...
                      gradient(gfr_vel_RA_act(:,3), 1/fs)];

    % Static Plots
    updateFigureContents('Position');
    grViz.plotPosition({estBody, actBody}, {'MIDPEL', 'LTIO', 'RTIO'});
    
    updateFigureContents('Animation Freeze');
    grid; view(0, 90); hold on;
    for i=idx(1):20:idx(end)
        grViz.plotLowerBody(estBody, i);
        grViz.plotLowerBody(actBody, i);
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
%         grViz.plotLowerBody(estBody, i);
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
%         grViz.plotLowerBody(actBody, i);
%         pause(1/1000);
%     end
% end