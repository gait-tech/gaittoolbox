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
    fnameV = 'TotalCaptureDataset/s1/acting1_BlenderZXY_YmZ.bvh';
    fnameS = 'TotalCaptureDataset/s1/Acting1_Xsens_AuxFields.sensors';
    
    % Check function input
    validateattributes(fnameV, {'char', 'struct'}, {});
    validateattributes(fnameS, {'char', 'struct'}, {});
    
    %% Initialization
    % Initialize variables and libraries
    addpath ./TCDlib;
    fs = 60;
    
    % Load calibration data
    calibIB = loadCalib('TotalCaptureDataset/IMUCalibration/calib_imu_bone.txt');
    calibIR = loadCalib('TotalCaptureDataset/IMUCalibration/calib_imu_ref.txt');
    
    % Load video orientation and position for each body segment
    if ischar(fnameV)
        dataV = loadBVHasAnimatedata(fnameV, 'mm');
    else
        dataV = fnameV;
    end     
    
    % Load sensor orientation and position for each body segment
    if ischar(fnameS)
        dataS = loadSensors(fnameS);
    else
        dataS = fnameS;
    end
    
    n = min(length(dataV.hips(1,:)), length(dataS.Pelvis(:,1)));
    keys = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    values = {'hips', 'lupleg', 'rupleg',...
        'lleg', 'rleg', 'lfoot', 'rfoot'};
    m = containers.Map(keys, values);
    for i=1:length(keys)
        dataS.(keys{i}) = dataS.(keys{i})(1:n,:);
        % dataV.(values{i}) = dataV.(values{i})(:,1:n)/1000;
        dataV.(values{i}) = quatrotate(calibIR.(keys{i}), dataV.(values{i})(:,1:n)'/1000)';
    end
    
    x0_pos_MP = dataV.hips(:,1)';
    x0_pos_LA = dataV.lfoot(:,1)';
    x0_pos_RA = dataV.rfoot(:,1)';
    x0_vel_MP = [0 0 0];
    x0_vel_LA = [0 0 0];
    x0_vel_RA = [0 0 0];
    sigma_acc = 0.5;
    
    gfr_acc_MP = quatrotate(quatconj(dataS.Pelvis(:,1:4)), dataS.Pelvis(:,5:7)) - [0 0 9.81];
    gfr_acc_LA = quatrotate(quatconj(dataS.L_LowLeg(:,1:4)), dataS.L_LowLeg(:,5:7)) - [0 0 9.81];
    gfr_acc_RA = quatrotate(quatconj(dataS.R_LowLeg(:,1:4)), dataS.R_LowLeg(:,5:7)) - [0 0 9.81];
    
    WIN_SECS = 0.25;
    VAR_WIN  = floor(fs*WIN_SECS); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;

    movVarAcc_pelvis = movingvar(sqrt( sum(dataS.Pelvis(:,5:7).^2,2)) ,VAR_WIN);
    bIsStatMP = movVarAcc_pelvis < 0;
    movVarAcc_lankle = movingvar(sqrt( sum(dataS.L_LowLeg(:,5:7).^2,2)) ,VAR_WIN);
    bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(dataS.R_LowLeg(:,5:7).^2,2)) ,VAR_WIN);
    bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
    
    %% -----------------------------------------------------------------------
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    uwb_mea.left_tibia_mid_pelvis = vecnorm((dataV.hips-dataV.lfoot)', 2, 2);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((dataV.hips-dataV.rfoot)', 2, 2);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.rfoot-dataV.lfoot)', 2, 2);

    [x_pri, x_pos] = kf_3_kmus(fs, sigma_acc, ...
        x0_pos_MP, x0_vel_MP, gfr_acc_MP, bIsStatMP,...
        x0_pos_LA, x0_vel_LA, gfr_acc_LA, bIsStatLA,...
        x0_pos_RA, x0_vel_RA, gfr_acc_RA, bIsStatRA, uwb_mea);
    
    % Visualise Trajectories for sanity check
    updateFigureContents('Animation Vicon Bone Segments');
    xlabel('x (m)');ylabel('y (m)');zlabel('z - Vertical (m)');
    hold on;grid on;axis('equal');view([-51 12]);%view([0 0]);
    
    % allIdx = 1:length(x_pos(:,1));
    allIdx = 1:300;
    
    % KF estimates
    bIsStationaryMP = repmat(bIsStatMP,1,3);
    bIsStationaryLA = repmat(bIsStatLA,1,3);
    bIsStationaryRA = repmat(bIsStatRA,1,3);
    [xMP,PxMP,yMP,PyMP,zMP,PzMP] = ...
        wrapper_PosVel_3KFs( fs,gfr_acc_MP,...
        [x0_pos_MP, x0_vel_MP]',1,bIsStationaryMP );

    [xLA,PxLA,yLA,PyLA,zLA,PzLA] = ...
        wrapper_PosVel_3KFs( fs,gfr_acc_LA,...
        [x0_pos_LA, x0_vel_LA]',1,bIsStationaryLA);

    [xRA,PxRA,yRA,PyRA,zRA,PzRA] = ...
        wrapper_PosVel_3KFs( fs,gfr_acc_RA,...
        [x0_pos_RA, x0_vel_RA]',1,bIsStationaryRA);
    
    plot3(xMP(allIdx,1),yMP(allIdx,1),zMP(allIdx,1),'.k');
    plot3(xLA(allIdx,1),yLA(allIdx,1),zLA(allIdx,1),'.b');
    plot3(xRA(allIdx,1),yRA(allIdx,1),zRA(allIdx,1),'.r');

    % KF with UWB
    plot3(x_pos(allIdx,1), x_pos(allIdx,2), x_pos(allIdx,3),'ok');
    plot3(x_pos(allIdx,7), x_pos(allIdx,8), x_pos(allIdx,9),'oc');
    plot3(x_pos(allIdx,13),x_pos(allIdx,14),x_pos(allIdx,15),'om');
    
    h_mp_mea = plot3(dataV.hips(1,allIdx),dataV.hips(2,allIdx),dataV.hips(3,allIdx),'-k');
    h_la_mea = plot3(dataV.lfoot(1,allIdx),dataV.lfoot(2,allIdx),dataV.lfoot(3,allIdx),'-b');
    h_ra_mea = plot3(dataV.rfoot(1,allIdx),dataV.rfoot(2,allIdx),dataV.rfoot(3,allIdx),'-r');
    
% end