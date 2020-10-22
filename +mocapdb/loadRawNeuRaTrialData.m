function data = loadRawNeuRaTrialData(dataDir, subjName, actName, ns, override, options)
    % Load data for selected movement trial
    % 
    % :param dataDir: data directory path
    % :param subjName: subject name of trial movement
    % :param actName: action name of trial movement
    % :param ns: yaw offset calibration mode. 
    %            NS1 = use yaw offset from ROM
    %            NS2 = use yaw offset from Vicon
    % :param override: [optional] override all calibration files
    % :param options: [optional] id to body mapping. By default uses the 
    %                 map specified in https://gait-tech.github.io/gaittoolbox/data.html.
    % 
    % :return: struct data with trial data inside
    %
    % .. Author: - Luke Wicent Sy (GSBME, Modified 2020 May 8)
    if nargin <= 4
        override = false;
    end
    if nargin <= 5
        options = struct('Pelvis', '00B40B91', ...
                        'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
                        'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
                        'L_LowLeg2', '00B40BA5', 'R_LowLeg2', '00B40C35', ...
                        'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
    end    
    ns = char(ns);
    DEGRANGE = (0:0.1:359) - 180;


    name = sprintf("%s-%s", subjName, actName);
    data = struct('name', name, ...
        'fnameV', sprintf('%s/rawvicon/%s-%s.mat', dataDir, subjName, actName), ...
        'fnameX', sprintf('%s/xsens/%s-%s.bvh', dataDir, subjName, actName), ...
        'fnameS', sprintf('%s/rawimu/%s-%s', dataDir, subjName, actName), ...
        'calibFnameSensorYawFixWorldFrame', ...
        sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrame.txt', ...
                dataDir, subjName, actName), ...
        'calibFnameSensorYawFixWorldFrameKFM', ...
        sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrameKFM.txt', ...
                dataDir, subjName, actName), ...
        'calibFnameSensorW2V', ...
        sprintf('%s/calib/%s-Calib-SensorW2V.txt', dataDir, subjName), ...
        'calibFnameSensorAccBias', ...
        sprintf('%s/calib/%s-%s-Calib-SensorAccBias.txt', dataDir, subjName, actName), ...
        'fnameRevStepDetect', ...
        sprintf('%s/rawstep-detect/%s-%s-revStepDetect.csv', ...
                dataDir, subjName, actName));

    data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);     
    % kinematic fitted model
    data.dataVKFM = mocapdb.ViconBody.loadViconMat(data.fnameV, false);
    data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
    data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
    data.dataS.fs = 100;
    data.calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
         sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dataDir, subjName), ...
         sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dataDir, subjName))' );
     
    % yaw off set calibration 
    if strcmp(ns(1:3), 'NS1')
        if exist(data.calibFnameSensorYawFixWorldFrame, 'file')
            data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        else
            % using ROM calibration
            dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/rawimu/%s-Trial-Walk-1', dataDir, subjName), options);
            dataS.fs = 100;           
            W__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/rawvicon/%s-Trial-Walk-1.mat', dataDir, subjName));
            W__dataV.changePosUnit('m', true);
            W__dataV = W__dataV.toWorldFrame(data.calibV2W);

            sIdx = max(W__dataV.getStartIndex()+1, 100);

            viconCalibSB = dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
            data.calibYawFix = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
            data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        end
    elseif strcmp(ns(1:3), 'NS2')
        if ~override && exist(data.calibFnameSensorYawFixWorldFrame, 'file')
            data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        else
            W__dataV = data.dataV.toWorldFrame(data.calibV2W);
            W__dataV.changePosUnit('m', true);
            data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
            data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        end
    elseif strcmp(ns(1:3), 'NS3')
        W__dataX = data.dataX;
        data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataX);
    else
        data.calibYawFix = struct();
        for j = ["Pelvis", "L_LowLeg", "R_LowLeg", "L_Foot", "R_Foot"]
            data.calibYawFix.(j).ori = [1 0 0 0];
        end
    end
    
    %% NS5
    if all(all(isnan(data.dataVKFM.PELV)))
        % all is nan
        data.calibYawFixKFM = [];
        fprintf('%s-%s KFM is NaN. Skipping generating SensorYawFixWorldFrame', subjName, actName);
    else
        if ~override && exist(data.calibFnameSensorYawFixWorldFrameKFM, 'file')
            data.calibYawFixKFM = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrameKFM);
        else
            W__dataV = data.dataVKFM.toWorldFrame(data.calibV2W);
            W__dataV.changePosUnit('m', true);
            data.calibYawFixKFM = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
            data.calibYawFixKFM.saveCalibCSV(data.calibFnameSensorYawFixWorldFrameKFM);
        end
    end
    
    if ~override && exist(data.calibFnameSensorW2V, 'file')
        data.calibW2V = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorW2V);
    else
          % using trackingmount calibration
%             data.calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
%                  sprintf('%s/calib/%s-Calib-SensorW2V.mat', dataDir, subjName), ...
%                  sprintf('%s/calib/%s-Calib-SensorW2V', dataDir, subjName), ...
%                  options, 100);
          % using ROM calibration
        dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/rawimu/%s-Trial-Walk-1', dataDir, subjName), options);
        dataS.fs = 100;           
        V__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/rawvicon/%s-Trial-Walk-1.mat', dataDir, subjName));
        V__dataV.changePosUnit('m', true);

        sIdx = max(V__dataV.getStartIndex()+1, 100);

        viconCalibSB = dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
        data.calibW2V = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
        data.calibW2V.saveCalibCSV(data.calibFnameSensorW2V);
    end

%     if false
%         % calculate acc pelvis bias from sitting trial.
%         % pelvis bias results to crouch divergence so not using it for now.
%         % not 100% tested but it works very well for the static trial
%         dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/rawimu/%s-Trial-TUG-1', dataDir, subjName), options);
%         dataS.fs = 100;               
%         dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/rawvicon/%s-Trial-TUG-1.mat', dataDir, subjName));
%         dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-TUG-1.bvh', dataDir, subjName), "mm");
%         data.bias = experiment.calcNeuRAPelvisAccBias01(dataS, dataV, ...
%                                 data.calibV2W, data.calibW2V, dataX, ...
%                                 100, 1000);
%     elseif false
%         % calculate acc pelvis bias from vicon data. Use static trial.
%         % pelvis bias results to crouch divergence so not using it for now.
%         % not 100% tested but it works very well for the static trial
%         dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/rawimu/%s-Trial-Static-1', dataDir, subjName), options);
%         dataS.fs = 100;               
%         dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/rawvicon/%s-Trial-Static-1.mat', dataDir, subjName));
%         dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-Static-1.bvh', dataDir, subjName), "mm");
%         data.bias = experiment.calcNeuRAPelvisAccBias02(dataS, dataV, ...
%                                 data.calibV2W, data.calibW2V, dataX, ...
%                                 1, -1);
%         data.bias
%     else
%    data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
%              'w__x', zeros(1, 3));
%     end
    data.revStepDetect = readtable(data.fnameRevStepDetect);
end

