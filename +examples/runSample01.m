 %% Simple example program 1
% Run ckf algorithm using `data/sample` comparing
%   1. ckf-3imu vs vicon system
%   2. xsens vs vicon system
%
% The sample data was an excerpt for the neura-sparse dataset

dir = 'data/sample';
expDir = sprintf('%s/output', dir);
mkdir(expDir);
mkdir(sprintf('%s/mat', dir));
addpath('mod-lib');

DEGRANGE = (0:0.1:359) - 180;
options = struct('Pelvis', '00B40B91', ...
                'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
                'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
                'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

ns = "NS2";
setups = {struct('est', 'ckf', 'accData', 'w__s', 'oriData', 'w__s', ...
               'initSrc', 'w__v', 'stepDetection', 'av03', ...
               'applyPred', 1, 'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
               'sigma2QAcc', 1e1^2); };
for i = 1:length(setups)
    setups{i}.label = getLabel(ns, setups{i});
end

n = struct('subj', 'S02', 'act', 'Trial-Walk-1', ...
           'startFrame', 1120, 'endFrame', 2074);

%% 1. ckf-3imu vs vicon system
name = sprintf("%s-%s-%s", ns, n.subj, n.act);
dataPath = sprintf('%s/mat/%s-%s-%s.mat', dir, ...
                    extractBetween(ns,1,3), n.subj, n.act);
if exist(dataPath, 'file')
    load(dataPath, 'data');
else
    data = struct('name', name, ...
        'fnameV', sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act), ...
        'fnameX', sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act), ...
        'fnameS', sprintf('%s/imu/%s-%s', dir, n.subj, n.act), ...
        'calibFnameSensorYawFixWorldFrame', ...
        sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrame.txt', ...
                dir, n.subj, n.act), ...
        'calibFnameSensorW2V', ...
        sprintf('%s/calib/%s-Calib-SensorW2V.txt', dir, n.subj), ...
        'fnameRevStepDetect', ...
        sprintf('%s/step-detect/%s-%s-revStepDetect.csv', ...
                dir, n.subj, n.act));

    data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);           
    if exist(data.fnameX, 'file')
        data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
    else
        data.dataX = [];
    end

    data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
    data.dataS.fs = 100;
    data.calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
         sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
         sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );

    if strcmp(extractBetween(ns,1,3), "NS2")
        if exist(data.calibFnameSensorYawFixWorldFrame, 'file')
            data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        else
            W__dataV = data.dataV.toWorldFrame(data.calibV2W);
            W__dataV.changePosUnit('m', true);
            data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
            data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        end
    else
        data.calibYawFix = struct();
        data.calibYawFix.Pelvis.ori = [1 0 0 0];
        data.calibYawFix.L_LowLeg.ori = [1 0 0 0];
        data.calibYawFix.R_LowLeg.ori = [1 0 0 0];
    end

    data.bias = struct('w__v', zeros(1, 3), 'w__x', zeros(1, 3));
    data.revStepDetect = readtable(data.fnameRevStepDetect);

    save(dataPath, 'data');
end

data.name = name;

fprintf("Data: %s\n", data.name);
r1 = examples.runNeuRASparse01Experiment(data.dataS, ...
        data.dataV, data.calibV2W, data.calibYawFix, ...
        data.dataX, data.revStepDetect, ...
        data.name, setups, expDir, n.startFrame, n.endFrame, data.bias);

%% 2. xsens vs vicon system
name = sprintf("%s-%s-%s", ns, n.subj, n.act);
load(sprintf("%s/%s-debug.mat", expDir, name));

sIdx = max(allIdx.w__v(1), allIdx.w__x(1));
eIdx = min(allIdx.w__v(end), allIdx.w__x(end));
nSamples = eIdx - sIdx + 1;

viconIdx0 = find(allIdx.w__v==sIdx,1):find(allIdx.w__v==eIdx,1);
xsensIdx0 = find(allIdx.w__x==sIdx,1):find(allIdx.w__x==eIdx,1);

csActBody = W__viconBody.getSubset(viconIdx0);
xsensBody = W__xsensBody.getSubset(xsensIdx0);

csActBodyRel = csActBody.changeRefFrame('MIDPEL');
estBodyRel = xsensBody.changeRefFrame('MIDPEL');
estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, xsensBody.qRPV);
csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);

results0a = xsensBody.diffRMSEandMean(csActBody);
results0 = estBody2.diffRMSEandMean(csActBody2);

results0.dPosW = results0a.dPos;
results0.name = name;
results0.label = sprintf("%s+viconvsxsens", ns);
results0.runtime = 0;
r2 = results0;
fprintf("%s/%s-%s\n", expDir, name, results0.label);
    

results = [struct2table(r1); struct2table(r2)];

function label = getLabel(ns, setup)    
    label = sprintf("%s+%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", ns, setup.est, ...
            setup.accData, setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyPred, ...
            setup.applyMeas, setup.applyCstr);
end