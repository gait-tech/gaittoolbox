%% Simple example program 2
% Run ckf algorithm using `data/sample` comparing
%   1. ckf-3imu vs vicon system
%   2. xsens vs vicon system
%
% The sample data was an excerpt for the neura-sparse dataset

%% Initialization
% set data source path
dir = 'data/sample';
% dir = 'data/neura-sparse01'; % uncomment when using neura-sparse dataset
% set result path
outDir = sprintf('%s/output', dir);
if ~exist(outDir)
    mkdir(outDir); 
end

dataList = { struct('subj', 'S02', 'act', 'Trial-Walk-1') };
% Uncomment when running the whole dataset 
% load data list from csv
% dataList = table2struct(readtable(sprintf('%s/data-list.csv', dir)));

% load additional libraries
addpath('mod-lib');

%% Loop through each trial listed in dataList
dataListN = length(dataList);
for i = 1:dataListN
    n = dataList{i};
    name = sprintf("%s-%s", n.subj, n.act);
    fprintf("Data %3d/%3d: %s\n", i, dataListN, data.name);
        
    %% Load trial data
    % vicon body in world frame
    dataV = mocapdb.ViconBody.loadCSV(sprintf("%s/vicon/%s", dir, name)); 
    % load imu measurements (acc, gyr, mag are in sensor frame, 
    % ori = orientation of sensor in world frame)
    [dataS, idx] = mocapdb.XsensBody.loadCSVs(sprintf("%s/imu/%s", dir, name)); 
    
    % change ori to orientation of body in world frame (w/ calibration)
    % note: not all IMUs have yaw fix offset
    % In the CKF 2020 paper, only the ankle yaw offset angles were used.
    calibYawFix = mocapdb.XsensBody.loadCalibCSV( ...
        sprintf('%s/calib/%s-Calib-SensorYawFixWorldFrame.txt', dir, name));
    calibYawFix.Pelvis.ori = [1 0 0 0];
    % apply yaw fix and obtain the sensor to body frame offset from the
    % initial vicon pose
    calibS2B = dataS.adjustFrame(calibYawFix, [1 0 0 0], true).calcCalibSB(dataV.togrBody(2, {}), 1);  
    % calibYawFix.qB * dataS.qB * calibS2B.qB for each body segment
    % dataB containts IMU measurement of body
    % dataB.Pelvis.ori = body in world frame
    % dataB.Pelvis.acc, gyr, mag = in body frame
    dataB = dataS.adjustFrame(calibYawFix, calibS2B.conj());
    % solve for gravity
    gfrAcc = dataB.calcGfrAcc();
    
    %% 1. ckf-3imu vs vicon system
    % generate reference grBody
    actBody = dataV.togrBody(1:dataV.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataV.fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
    % load step detect data
    revStepDetect = readtable(sprintf('%s/step-detect/%s-revStepDetect.csv', ...
                    dir, name));
    revStepDetect = revStepDetect(idx(1):idx(2),:);
                
    % calculate minimum knee angle from actual data.
    % since the person is assumed to be standing straight at start, we set
    % minimum angle to starting knee angle OR zero (if negative)
    alphaLKmin = actBody.calcJointAnglesLKnee(1);
    alphaLKmin = min(alphaLKmin(2), 0);
    alphaRKmin = actBody.calcJointAnglesRKnee(1);
    alphaRKmin = min(alphaRKmin(2), 0);

    % calculate body length
    body = struct('PV_d', actBody.calcPelvisLength(), ...
                  'LT_d', actBody.calcLFemurLength(), ...
                  'RT_d', actBody.calcRFemurLength(), ...
                  'LS_d', actBody.calcLShankLength(), ...
                  'RS_d', actBody.calcRShankLength() );

    % get initial state from vicon grBody
    vel = actBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
    x0 = [actBody.MIDPEL(1,:) vel.MIDPEL(1,:) ...
          actBody.LTIO(1,:) vel.LTIO(1,:) ...
          actBody.RTIO(1,:) vel.RTIO(1,:) ]';  

    dt = 1.0/dataS.fs;
    x0(4:6, :) = (x0(4:6, :)*dt - 0.5*gfrAcc.Pelvis(1,:)'*dt^2)/dt;
    x0(10:12, :) = (x0(10:12, :)*dt - 0.5*gfrAcc.L_LowLeg(1,:)'*dt^2)/dt;
    x0(16:18, :) = (x0(16:18, :)*dt - 0.5*gfrAcc.R_LowLeg(1,:)'*dt^2)/dt;

    v3Options = struct('fs', dataS.fs, 'applyMeas', 76, 'applyCstr', 355, ...
                       'sigma2QAccMP', 1e1^2, 'sigma2QAccLA', 1e1^2, ...
                       'sigma2QAccRA', 1e1^2, ...
                       'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);
                
    t0 = cputime;
    [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.ckf_3imus( x0, 0.5, ...
                    gfrAcc.Pelvis, false(actBody.nSamples), dataB.Pelvis.ori, ...
                    gfrAcc.L_LowLeg, revStepDetect.stepL, dataB.L_LowLeg.ori, ...
                    gfrAcc.R_LowLeg, revStepDetect.stepR, dataB.R_LowLeg.ori, ...
                    body.PV_d, body.LT_d, body.RT_d, body.LS_d, body.RS_d, ...
                    v3Options);
    runtime = cputime-t0;
    
    % generate grBody from ckf filter output
    estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
                   'xyzColor', {'r', 'g', 'b'}, 'fs', dataS.fs, ...
                   'MIDPEL', x_pos_v2(:, 1:3), ...
                   'LFEP', t_dat_v2.LFEP, ...
                   'LFEO', t_dat_v2.LFEO, ...
                   'LTIO', x_pos_v2(:, 7:9), ...
                   'RFEP', t_dat_v2.RFEP, ...
                   'RFEO', t_dat_v2.RFEO, ...
                   'RTIO', x_pos_v2(:, 13:15), ...
                   'qRPV', dataB.Pelvis.ori, ...
                   'qLTH', t_dat_v2.qLTH(:, :), ...
                   'qRTH', t_dat_v2.qRTH(:, :), ...
                   'qLSK', dataB.L_LowLeg.ori, ...
                   'qRSK', dataB.R_LowLeg.ori );
               
    % convert grBody to midpel frame
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    % revert back to global frame but with MIDPEL position equal to actBody
    estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, estBody.qRPV);
 
    % compare grBody of vicon (actBody) and estimate (estBody)
    r1 = estBody2.diffRMSEandMean(actBody);
    r1.name = name;
    r1.label = "ckf";
    r1.runtime = 0;
    
    %% 2. xsens vs vicon system
    dataX = mocapdb.BVHBody.loadXsensBVHFile( ...
            sprintf('%s/xsens/%s.bvh', dir, name), "mm");
    % transformation to express dataX in biomechanical convention
    % x = forward, z = upward, y = left
    qXsens2BMC = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
    xsnBody = dataX.getSubset(idx(1):idx(2)).toWorldFrame(qXsens2BMC) ...
                .changePosUnit('m') ...
                .togrBody(1:dataS.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataS.fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
    
    % convert grBody to midpel frame
    xsnBodyRel = xsnBody.changeRefFrame('MIDPEL');
    % revert back to global frame but with MIDPEL position equal to actBody
    xsnBody2 = xsnBodyRel.toWorldFrame(actBody.MIDPEL, xsnBody.qRPV);
 
    % compare grBody of vicon (actBody) and xsens (xsnBody)
    r2 = xsnBody2.diffRMSEandMean(actBody);
    r2.name = name;
    r2.label = "viconvsxsens";
    r2.runtime = 0;
end    

results = [struct2table(r1); struct2table(r2)];