% ======================================================================
%> Generate step detect templates from variance
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
outDir = sprintf('%s/step-detect', dir);
ns = "NS2";

dataList = readtable(sprintf('%s/data-list-v2.csv', dir));
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
stepDetectWindow = 0.25;
stepDetectThreshold = 1;

results = table();
dataN = size(dataList, 1);

for i = 1:dataN
    n = table2struct(dataList(i, :));

    name = sprintf("%s-%s", n.subj, n.act);
    fnameV = sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act);
    fnameX = sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act);
    fnameS = sprintf('%s/imu/%s-%s', dir, n.subj, n.act);
    calibFnameSensorYawFixWorldFrame = sprintf('%s/calib/%s-Calib-SensorYawFixWorldFrame.txt', dir, n.subj);
    calibFnameSensorW2V = sprintf('%s/calib/%s-Calib-SensorW2V.txt', dir, n.subj);
    
    fs = 100;
    dataV = mocapdb.ViconBody.loadViconMat(fnameV);
    dataX = mocapdb.BVHBody.loadXsensBVHFile(fnameX, "mm");
    dataS = mocapdb.XsensBody.loadMTExport(fnameS, options);
    dataS.fs = 100;
    calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
             sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
             sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );
    
    W__dataV = dataV.toWorldFrame(calibV2W);
    W__dataV.changePosUnit('m', true);
    calibYawFix = dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);

    bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                  'w__x', zeros(1, 3));

    display(sprintf("Data %3d/%3d: %s", i, dataN, name));

    W__dataV = dataV.toWorldFrame(calibV2W);
    W__dataV.changePosUnit('m', true);
    W__dataS = dataS;
    W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__v;
    % apply yaw offset to orientation
    W__dataS.L_LowLeg.ori = quatmultiply(calibYawFix.L_LowLeg.ori, W__dataS.L_LowLeg.ori);
    W__dataS.R_LowLeg.ori = quatmultiply(calibYawFix.R_LowLeg.ori, W__dataS.R_LowLeg.ori);

    %% generate step detection file from Vicon acceleration
    W__viconBody = W__dataV.togrBody(1:dataV.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                     'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                     'xyzColor', {'m', 'y', 'c'}});  
    vel = W__viconBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
    acc = W__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});   
    
    VAR_WIN  = floor(fs*stepDetectWindow); % NUM_SAMPLES
    ACC_VAR_THRESH = stepDetectThreshold;

    movVarAcc_lankle = movingvar(sqrt( sum(acc.LTIO .^2, 2)), VAR_WIN);
    stepL = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(acc.RTIO .^2, 2)), VAR_WIN);
    stepR = movVarAcc_rankle < ACC_VAR_THRESH;
        
    stepDetectVicon = table(stepL, stepR);
    writetable(stepDetectVicon, sprintf('%s/%s-viconStepDetect.csv', outDir, name));
    
    %% generate step detection file from IMU measured acceleration
    accS = {};
    accS.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                            W__dataS.Pelvis.acc) - [0 0 9.81];
    accS.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                            W__dataS.L_LowLeg.acc) - [0 0 9.81];
    accS.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                            W__dataS.R_LowLeg.acc) - [0 0 9.81];

    VAR_WIN  = floor(fs*stepDetectWindow); % NUM_SAMPLES
    ACC_VAR_THRESH = stepDetectThreshold;

    movVarAcc_lankle = movingvar(sqrt( sum(accS.LA .^2, 2)), VAR_WIN);
    stepL = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(accS.RA .^2, 2)), VAR_WIN);
    stepR = movVarAcc_rankle < ACC_VAR_THRESH;
    
    stepDetectIMU = table(stepL, stepR);
    writetable(stepDetectIMU, sprintf('%s/%s-imuStepDetect.csv', outDir, name));
end