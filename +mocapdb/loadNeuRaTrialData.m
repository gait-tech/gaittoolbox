function data = loadNeuRaTrialData(dataDir, subjName, actName, usebuffer, prefix, removeaccbias)
    % Load data for selected movement trial
    % 
    % :param dataDir: data directory path
    % :type dataDir: characters, string
    % :param subjName: subject name of trial movement
    % :type subjName: characters, string
    % :param actName: action name of trial movement
    % :type actName: characters, string
    % :param usebuffer: use buffer .mat if available
    % :type usebuffer: Optional, boolean
    % :param prefix: prefix to name
    % :type prefix: Optional, characters, string
    % 
    % :return: struct data with trial data inside
    % :rtype: struct
    %
    % .. Author: - Luke Wicent Sy (GSBME, Modified 2020 May 25)
    
    if nargin <= 3
        usebuffer = true;
    end
    if nargin <= 4
        prefix = '';
    end
    
    name = sprintf("%s%s-%s", prefix, subjName, actName);
    dataPath = sprintf('%s/mat/%s%s-%s-buffer.mat', dataDir, ...
                        prefix, subjName, actName);
    if usebuffer && exist(dataPath, 'file')
        load(dataPath, 'data');
    else
        data = struct('name', name, ...
            'fnameV', sprintf('%s/vicon/%s-%s', dataDir, subjName, actName), ...
            'fnameVKFM', sprintf('%s/vicon-kfm/%s-%s', dataDir, subjName, actName), ...
            'fnameX', sprintf('%s/xsens/%s-%s.bvh', dataDir, subjName, actName), ...
            'fnameS', sprintf('%s/imu/%s-%s', dataDir, subjName, actName), ...
            'calibFnameSensorYawFixWorldFrame', ...
            sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrame.txt', ...
                    dataDir, subjName, actName), ...
            'calibFnameSensorYawFixWorldFrameKFM', ...
            sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrameKFM.txt', ...
                    dataDir, subjName, actName), ...
            'calibFnameSensorAccBias', ...
            sprintf('%s/calib/%s-%s-Calib-SensorAccBias.txt', dataDir, subjName, actName), ...
            'calibFnameSensorW2V', ...
            sprintf('%s/calib/%s-Calib-SensorW2V.txt', dataDir, subjName), ...
            'fnameRevStepDetect', ...
            sprintf('%s/step-detect/%s-%s-revStepDetect.csv', ...
                    dataDir, subjName, actName));

        data.dataV = mocapdb.ViconBody.loadCSV(data.fnameV); 
        data.dataVKFM = mocapdb.ViconBody.loadCSV(data.fnameVKFM); 
        % load imu measurements (acc, gyr, mag are in sensor frame, 
        % ori = orientation of sensor in world frame)
        [data.dataS, data.idx] = mocapdb.XsensBody.loadCSVs(data.fnameS); 
        data.dataV.ftStartIndex = data.idx(1);
        data.dataV.ftEndIndex = data.idx(2);
        
        data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
        % transformation to express dataX in biomechanical convention
        % x = forward, z = upward, y = left
        qXsens2BMC = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        data.dataX = data.dataX.getSubset(data.idx(1):data.idx(2)) ...
                        .toWorldFrame(qXsens2BMC).adjustFootFrame2BMC() ...
                        .changePosUnit('m');        
    
        data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        if exist(data.calibFnameSensorYawFixWorldFrameKFM, 'file')
            data.calibYawFixKFM = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrameKFM);
        else
            data.calibYawFixKFM = [];
        end
            
        data.calibW2V = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorW2V);
        data.calibS2B = [];
        % data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
        %                    'w__x', zeros(1, 3));
        data.bias = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorAccBias);
        data.revStepDetect = readtable(data.fnameRevStepDetect);
        
        if exist(dataDir, 'dir')
            if ~exist(sprintf('%s/mat', dataDir), 'dir')
                mkdir(sprintf('%s/mat', dataDir));
            end
            save(dataPath, 'data');
        end
    end
end

