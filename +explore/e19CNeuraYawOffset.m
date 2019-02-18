global DEGRANGE;
DEGRANGE = (0:0.1:359) - 180;

dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);

dataList = readtable(sprintf('%s/data-list-v2.csv', dir));
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

dataN = size(dataList, 1);
results = table();

for i = 1:dataN
    n = table2struct(dataList(i, :));
    
    name = sprintf("%s-%s-%s", 'neura', n.subj, n.act);
    dataPath = sprintf('%s/mat/%s.mat', dir, name);

    calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
             sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
             sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );
         
    calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
                 sprintf('%s/calib/%s-Calib-SensorW2V.mat', dir, n.subj), ...
                 sprintf('%s/calib/%s-Calib-SensorW2V', dir, n.subj), ...
                 options, 100); 
             
    calibYawFix = mocapdb.XsensBody.loadCalibCSV( ...
                 sprintf('%s/calib/%s-Calib-SensorYawFixWorldFrame.txt', dir, n.subj));
%     rad2deg(quat2eul(calibW2V.Pelvis.ori))

    dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-%s', dir, n.subj, n.act), options);
    dataS.fs = 100;
    dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act));
    fs = dataV.fs;

    nSamples = min(dataV.nSamples, dataS.nSamples);
    V__dataV = dataV.getSubset(1:nSamples);
    V__dataV.changePosUnit('m', true);
    W__dataS = dataS.getSubset(1:nSamples);
    % V__dataS = W__dataS.toViconFrame(calibW2V);
    
    V__dataV = V__dataV.toWorldFrame(calibV2W);
    V__dataS = W__dataS;

    sIdx = max(V__dataV.getStartIndex()+1, 100);
    eIdx = length(V__dataV.PELV(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);

    V__viconBody = V__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', V__dataV.fs, ...
                             'xyzColor', {'m', 'y', 'c'}});  
    refAcc = V__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
    refAcc.MIDPEL = refAcc.MIDPEL(sIdx:eIdx,:);
    refAcc.LTIO = refAcc.LTIO(sIdx:eIdx,:);
    refAcc.RTIO = refAcc.RTIO(sIdx:eIdx,:);

    sensors = dataS.exportRawMeasurementAsStruct({'Pelvis', 'L_LowLeg', 'R_LowLeg'}, ...
                        {'PELV', 'LANK', 'RANK'});

    estAcc = {};
    estAcc.MIDPEL = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                            W__dataS.Pelvis.acc) - [0 0 9.81];
    estAcc.MIDPEL = estAcc.MIDPEL(sIdx:eIdx,:);
    % estAcc.MIDPEL = quatrotate(quatconj(calibW2V.Pelvis.ori), estAcc.MIDPEL);
    estAcc.LTIO = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                            W__dataS.L_LowLeg.acc) - [0 0 9.81];
    estAcc.LTIO = estAcc.LTIO(sIdx:eIdx,:);
    % estAcc.LTIO = quatrotate(quatconj(calibW2V.L_LowLeg.ori), estAcc.LTIO);
    estAcc.RTIO = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                            W__dataS.R_LowLeg.acc) - [0 0 9.81];
    estAcc.RTIO = estAcc.RTIO(sIdx:eIdx,:);
    % estAcc.RTIO = quatrotate(quatconj(calibW2V.R_LowLeg.ori), estAcc.RTIO);

    [thetaP, errP] = findOptimalThetaBrute(refAcc.MIDPEL, estAcc.MIDPEL);
    [thetaL, errL] = findOptimalThetaBrute(refAcc.LTIO, estAcc.LTIO);
    [thetaR, errR] = findOptimalThetaBrute(refAcc.RTIO, estAcc.RTIO);

    [~, peakIdxL] = min(abs(DEGRANGE-thetaL)); 
    [~, peakIdxR] = min(abs(DEGRANGE-thetaR)); 
    
    calibThetaL = rad2deg(quat2eul(quatconj(calibYawFix.L_LowLeg.ori)));
    calibThetaR = rad2deg(quat2eul(quatconj(calibYawFix.R_LowLeg.ori)));

    
    results0 = {name thetaP thetaL thetaR calibThetaL(1) calibThetaR(1)};
    results = [results; results0];
    display(sprintf("%3d/%3d", i, dataN));
end

results.Properties.VariableNames = {'name', 'viconP', 'viconL', 'viconR', 'ROML', 'ROMR'};
writetable(results, 'neura-sparse01/explore-v2/yawoffsets20190207.csv')

function [theta, err] = findOptimalThetaBrute(refAcc, estAcc)
    global DEGRANGE;
    err = zeros(length(DEGRANGE), 1);
    for i=1:length(DEGRANGE)
        estAcc2 = rotateThetaAlongZ(estAcc, DEGRANGE(i));
        err(i) = sqrt(mean(nanmean((refAcc(:,1:2) - estAcc2(:,1:2)).^2, 1)));
    end
    [M, I] = min(err);
    theta = DEGRANGE(I);
end

function newData = rotateThetaAlongZ(data, deg)
    R = [cosd(deg) sind(deg) 0;
         -sind(deg) cosd(deg) 0;
             0 0 1];
    newData = quatrotate(quatconj(rotm2quat(R)), data);
end