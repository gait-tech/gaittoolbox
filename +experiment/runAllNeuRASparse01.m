% ======================================================================
%> Run experiment for all NeuRA data
%> NS1: use yaw offset from ROM
%> NS2: use yaw offset from Vicon
%> writetable(results, 'neura-sparse01/explore-v2/results20190204.csv')
% ======================================================================
dir = 'neura-sparse01';
% expDir = sprintf('%s/explore', dir);
expDir = sprintf('%s/explore-v2', dir);

DEGRANGE = (0:0.1:359) - 180;
dataList = readtable(sprintf('%s/data-list-v2.csv', dir));

% options = struct('Pelvis', '00B40B91', ...
%     'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%     'L_LowLeg', '00B40C49', 'R_LowLeg', '00B40C4A', ...
%     'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
%     'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
% options = struct('Pelvis', 'Pelvis', ...
%     'L_UpLeg', 'LeftUpperLeg', 'R_UpLeg', 'RightUpperLeg', ...
%     'L_LowLeg', 'prop', 'R_LowLeg', 'prop_1', ...
%     'L_Foot', 'LeftFoot', 'R_Foot', 'RightFoot');

results = table();

nsList = {'NS2'};
% nsIdx = 1;
% for i=[0.0:0.01:0.1 0.12:0.02:0.3]
%     nsList{nsIdx} = sprintf('NS2d%0.2f', i);
%     nsIdx = nsIdx + 1;
% end

for nsI = 1:length(nsList)
    ns = nsList{nsI};
    
    setups = {
    };
    setups{end+1} = struct('est', 'ekfv3', ...
               'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
               'initSrc', 'w__v', 'stepDetection', 'av03', ...
               'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
               'sigmaQAcc', 1e1);
%     if strcmp(ns, 'NS2')
%         mIList = [302, 376];
%         cIList = [0 1 351 355];
%         setups{end+1} = struct('est', 'ekfv3', ...
%                            'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
%                            'initSrc', 'w__v', 'stepDetection', 'av03', ...
%                            'applyMeas', 2, 'applyCstr', 355, 'P', 0.5, ...
%                            'sigmaQAcc', 1e1);
%         setups{end+1} = struct('est', 'ekfv3', ...
%                            'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
%                            'initSrc', 'w__v', 'stepDetection', 'av03', ...
%                            'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
%                            'sigmaQAcc', 1e1);
%     else
%         mIList = 302;
%         cIList = 351;
%         
%         setups{end+1} = struct('est', 'ekfv3', ...
%                            'accData', 'w__v', 'oriData', 'w__v', 'accDataNoise', 0, ...
%                            'initSrc', 'w__v', 'stepDetection', 'av03', ...
%                            'applyMeas', 0, 'applyCstr', 2, 'P', 0.5, ...
%                            'sigmaQAcc', 1e1, ...
%                            'sigmaUwbLeg', 0);
%     end

    for pI = [1 11]
        for mI = [11]
            for cI = [7]
                for sdI = {'av03'} % {'av01', 'av03'}
    %                     setups{end+1} = struct('est', 'ekfv3', ...
    %                            'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
    %                            'initSrc', 'w__v', 'stepDetection', sdI, ...
    %                            'applyMeas', mI, 'applyCstr', cI, 'P', 0.5, ...
    %                            'sigmaQAcc', 1e1, ...
    %                            'sigmaUwbLeg', suwb);
                    setups{end+1} = struct('est', 'lieekfv1', ...
                               'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                               'initSrc', 'w__v', 'stepDetection', sdI, ...
                               'applyPred', pI, 'applyMeas', mI, ...
                               'applyCstr', cI, 'P', 0.5, ...
                               'sigmaQAcc', 1e1, 'sigmaQAngVel', 1e2);
                end
            end
        end
    end

    for i = 1:length(setups)
        setups{i}.label = getLabel(ns, setups{i});
    end

    dataN = size(dataList, 1);

    for i = [15]
        n = table2struct(dataList(i, :));
        
        uwbDistSigma = 0.0;
        if strcmp(ns(1:3), 'NS2') && (size(ns, 2) > 3)
            uwbDistSigma = str2double(ns(5:end)); 
        end
        
        if ~(strcmp(n.act(7:end-2), "SpeedSkater") || ...
           strcmp(n.act(7:end-2), "Static") || ...
           strcmp(n.act(7:end-2), "HighKneeJog")) && ...
           uwbDistSigma <= (0.03+1e-4)
%             setups{1}.sigmaUwbLeg = 1e-1;
            setups{1}.sigmaUwbLeg = 1e0;
        else
            setups{1}.sigmaUwbLeg = 1e0;
            continue;
        end

        name = sprintf("%s-%s-%s", ns, n.subj, n.act);
        dataPath = sprintf('%s/mat/%s-%s-%s.mat', dir, ns(1:3), n.subj, n.act);
        if exist(dataPath, 'file')
            load(dataPath, 'data');
        else
            data = struct('name', name, ...
                'fnameV', sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act), ...
                'fnameX', sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act), ...
                'fnameS', sprintf('%s/imu/%s-%s', dir, n.subj, n.act), ...
                'calibFnameSensorYawFixWorldFrame', ...
                sprintf('%s/calib/%s-Calib-SensorYawFixWorldFrame.txt', ...
                        dir, n.subj), ...
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
            
            if strcmp(n.subj, 'S03')
                nSamples = min(data.dataV.nSamples, data.dataX.nSamples);
                data.dataV = data.dataV.getSubset(1:nSamples);
                data.dataX = data.dataX.getSubset(1:nSamples);
                
                qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
                dataXtmp = data.dataX.toWorldFrame(qXsensV2W);
                data.dataV.qLTH = dataXtmp.qLeftUpLeg;
                data.dataV.qLSK = dataXtmp.qLeftLeg;
%                 sIdx = max(1, n.startFrame);
%                 qOffset = quatmultiply(dataXtmp.qLeftLeg(sIdx,:), ...
%                                 quatconj(data.dataV.qLSK(sIdx,:)));
%                 data.dataV.qLSK = quatmultiply(qOffset, data.dataV.qLSK);
%                 data.dataV.qLTH = quatmultiply(axang2quat([0 0 1 deg2rad(-50)]), data.dataV.qLTH);
%                 data.dataV.qLSK = quatmultiply(axang2quat([0 0 1 deg2rad(-50)]), data.dataV.qLSK);
                dLTH = vecnorm(data.dataV.LFEO-data.dataV.LFEP, 2, 2);
                dLSK = vecnorm(data.dataV.LFEO-data.dataV.LTIO, 2, 2);
                LTHz = quat2rotm(data.dataV.qLTH); 
                LTHz = squeeze(LTHz(:,3,:))';
                LSKz = quat2rotm(data.dataV.qLSK); 
                LSKz = squeeze(LSKz(:,3,:))';
                data.dataV.LFEO = data.dataV.LFEP - dLTH.*LTHz;
                data.dataV.LTIO = data.dataV.LFEO - dLSK.*LSKz;
            end
            
            data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
            data.dataS.fs = 100;
            data.calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
                 sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
                 sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );
            if strcmp(ns(1:3), 'NS1')
                if exist(data.calibFnameSensorYawFixWorldFrame, 'file')
                    data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                else
                    % using ROM calibration
                    dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Walk-1', dir, n.subj), options);
                    dataS.fs = 100;           
                    W__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Walk-1.mat', dir, n.subj));
                    W__dataV.changePosUnit('m', true);
                    W__dataV = W__dataV.toWorldFrame(data.calibV2W);

                    sIdx = max(W__dataV.getStartIndex()+1, 100);

                    viconCalibSB = dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
                    data.calibYawFix = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
                    data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                end
            else %strcmp(ns(1:3), 'NS2')
                W__dataV = data.dataV.toWorldFrame(data.calibV2W);
                W__dataV.changePosUnit('m', true);
                data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
            end

            if exist(data.calibFnameSensorW2V, 'file')
                data.calibW2V = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorW2V);
            else
                  % using trackingmount calibration
    %             data.calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
    %                  sprintf('%s/calib/%s-Calib-SensorW2V.mat', dir, n.subj), ...
    %                  sprintf('%s/calib/%s-Calib-SensorW2V', dir, n.subj), ...
    %                  options, 100);
                  % using ROM calibration
                dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Walk-1', dir, n.subj), options);
                dataS.fs = 100;           
                V__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Walk-1.mat', dir, n.subj));
                V__dataV.changePosUnit('m', true);

                sIdx = max(V__dataV.getStartIndex()+1, 100);

                viconCalibSB = dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
                data.calibW2V = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
                data.calibW2V.saveCalibCSV(data.calibFnameSensorW2V);
            end

%             if false
%                 % calculate acc pelvis bias from sitting trial.
%                 % pelvis bias results to crouch divergence so not using it for now.
%                 % not 100% tested but it works very well for the static trial
%                 dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-TUG-1', dir, n.subj), options);
%                 dataS.fs = 100;               
%                 dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-TUG-1.mat', dir, n.subj));
%                 dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-TUG-1.bvh', dir, n.subj), "mm");
%                 data.bias = experiment.calcNeuRAPelvisAccBias01(dataS, dataV, ...
%                                         data.calibV2W, data.calibW2V, dataX, ...
%                                         100, 1000);
%             elseif false
%                 % calculate acc pelvis bias from vicon data. Use static trial.
%                 % pelvis bias results to crouch divergence so not using it for now.
%                 % not 100% tested but it works very well for the static trial
%                 dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Static-1', dir, n.subj), options);
%                 dataS.fs = 100;               
%                 dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Static-1.mat', dir, n.subj));
%                 dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-Static-1.bvh', dir, n.subj), "mm");
%                 data.bias = experiment.calcNeuRAPelvisAccBias02(dataS, dataV, ...
%                                         data.calibV2W, data.calibW2V, dataX, ...
%                                         1, -1);
%                 data.bias
%             else
                data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                          'w__x', zeros(1, 3));
%             end
            
            data.revStepDetect = readtable(data.fnameRevStepDetect);
            
            save(dataPath, 'data');
        end

        uwbDistSigma = 0.0;
        if strcmp(ns(1:3), 'NS2') && (size(ns, 2) > 3)
            uwbDistSigma = str2double(ns(5:end)); 
        end
        data.name = name;
        
        fprintf("Data %3d/%3d: %s\n", i, dataN, data.name);
        r = experiment.runNeuRASparse01Experiment(data.dataS, ...
                data.dataV, data.calibV2W, data.calibYawFix, data.calibW2V, ...
                data.dataX, data.revStepDetect, uwbDistSigma, ...
                data.name, setups, expDir, n.startFrame, n.endFrame, data.bias);
        results = [results; struct2table(r)];
    end
end

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
    if setup.accData == 'v'
        if setup.accDataNoise == 0 
            aD = 'v';
        else
            aD = strrep(sprintf('v%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData;
    end
    if strcmp(setup.est, 'ekfv3')
        label = sprintf("%s+A%sO%sI%s+S%s+M%02d+C%03d", ns, aD, ...
            setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyMeas, setup.applyCstr);
    else
        label = sprintf("%s+%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", ns, setup.est, ...
            aD, setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyPred, ...
            setup.applyMeas, setup.applyCstr);
    end
end