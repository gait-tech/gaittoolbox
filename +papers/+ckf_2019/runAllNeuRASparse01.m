% Run experiment for all NeuRA data
% NS2: use yaw offset from Vicon
%
% Example:
%       writetable(results, 'neura-sparse01/explore/results20190812.csv')
%

dir = 'data/neura-sparse01';
expDir = sprintf('%s/output', dir);
mkdir(expDir);
mkdir(sprintf('%s/mat', dir));
addpath('mod-lib');

DEGRANGE = (0:0.1:359) - 180;
dataList = readtable('+papers/+ckf_2019/data-list.csv');
dataN = size(dataList, 1);
DATARANGE = 1:dataN;

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

results = table();

for ns = ["NS2"]   
    setups = {struct('est', 'ckf', ...
                   'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                   'initSrc', 'w__v', 'stepDetection', 'av03', ...
                   'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
                   'sigma2QAcc', 1e1^2); };

    for i = 1:length(setups)
        setups{i}.label = getLabel(ns, setups{i});
    end

    for i = DATARANGE
        n = table2struct(dataList(i, :));
        
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
        
        fprintf("Data %3d/%3d: %s\n", i, dataN, data.name);
        r = papers.ckf_2019.runNeuRASparse01Experiment(data.dataS, ...
                data.dataV, data.calibV2W, data.calibYawFix, ...
                data.dataX, data.revStepDetect, ...
                data.name, setups, expDir, n.startFrame, n.endFrame, data.bias);
        results = [results; struct2table(r)];
    end
end

rIdx = size(results, 1) + 1;
results = table2struct(results);

%% file list vicon vs xsens comparison
for i = DATARANGE
    n = table2struct(dataList(i, :));
    name = sprintf("%s-%s-%s", ns, n.subj, n.act);
    load(sprintf("%s/%s-debug.mat", expDir, name));
    
    sIdx = max(allIdx.w__v(1), allIdx.w__x(1));
    eIdx = min(allIdx.w__v(end), allIdx.w__x(end));
    nSamples = eIdx - sIdx + 1;
    
    viconIdx0 = find(allIdx.w__v==sIdx,1):find(allIdx.w__v==eIdx,1);
    xsensIdx0 = find(allIdx.w__x==sIdx,1):find(allIdx.w__x==eIdx,1);
    
    csActBody = W__viconBody.getSubset(viconIdx0);
    estBody = W__xsensBody.getSubset(xsensIdx0);
       
    csActBodyRel = csActBody.changeRefFrame('MIDPEL');
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
    csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    
    results0a = estBody.diffRMSEandMean(csActBody);
    results0 = estBody2.diffRMSEandMean(csActBody2);
    
    results0.dPosW = results0a.dPos;
    results0.name = name;
    results0.label = sprintf("%s+viconvsxsens", ns);
    results0.runtime = 0;
    results(rIdx) = results0;
    rIdx = rIdx + 1;
    fprintf("%s/%s-%s\n", expDir, name, results0.label);
    
%     targetname = sprintf('%s/%s-viconvsxsens', outDir, name);
%     estBody.exportc3d(sprintf('%s.c3d', targetname), struct(), csActBody);
end

results = struct2table(results);

% Append new results
dataPath = sprintf("%s/results.mat", expDir);
if exist(dataPath, 'file')
    newResults = results;
    load(dataPath);
    [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
    results(ia,:) = [];
    results = [results; newResults];
end
save(sprintf("%s/results.mat", expDir), 'results')

function label = getLabel(ns, setup)    
    label = sprintf("%s+%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", ns, setup.est, ...
            setup.accData, setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyPred, ...
            setup.applyMeas, setup.applyCstr);
end