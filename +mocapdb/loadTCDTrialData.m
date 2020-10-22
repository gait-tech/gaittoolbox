function data = loadTCDTrialData(dataDir, subjName, actName, usebuffer, prefix)
    % Load data for selected movement trial
    % 
    % :param dataDir: data directory path
    % :type dataDir: characters, string
    % :param subjName: subject name of trial movement
    % :type subjName: characters, string
    % :param actName: action name of trial movement
    % :type actName: characters, string
    % :param usebuffer: use buffer .mat if available
    % :type usebuffer: boolean, optional
    % :param prefix: prefix to name
    % :type prefix: characters, string
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
        actNameL = lower(actName);
        data = struct('name', name, ...
            'fnameV', sprintf('%s/vicon/%s/%s_BlenderZXY_YmZ.bvh', dataDir, subjName, actNameL), ...
            'fnameS', sprintf('%s/gyroMag/%s/%s_Xsens_AuxFields.sensors', dataDir, subjName, actName), ...
            'calibFnameSensor2Body', sprintf('%s/imu/%s/%s_%s_calib_imu_bone.txt', dataDir, subjName, subjName, actNameL), ...
            'calibFnameSensorW2V', sprintf('%s/imu/%s/%s_%s_calib_imu_ref.txt', dataDir, subjName, subjName, actNameL), ...
            'fnameRevStepDetect', sprintf('%s/step-detect/%s-%s-revStepDetect.csv', dataDir, subjName, actName) ...
            );

        % load vicon (optical motion tracker) data
        data.calibV2W = rotm2quat([1 0 0; 0 0 -1; 0 1 0]);
        data.dataV = mocapdb.BVHBody.loadBVHFile(data.fnameV, 'mm');
        data.dataV = data.dataV.toWorldFrame(data.calibV2W) ...
                        .adjustFootFrame2BMC().toViconBody();
        data.dataV.changePosUnit('m', true);
        data.dataV.fs = 60;
        
        % load imu measurements (acc, gyr, mag are in sensor frame, 
        % ori = orientation of sensor in world frame)
        data.dataS = mocapdb.XsensBody.loadSensorFile(data.fnameS);
        data.dataS.fs = 60;
        
        data.nSamples = min(data.dataV.nSamples, data.dataS.nSamples);
        if data.dataV.nSamples ~= data.dataS.nSamples
            fprintf("dataV (%d) and dataS (%d) have different nSamples. Set nSamples=%d\n", ...
                    data.dataV.nSamples, data.dataS.nSamples, data.nSamples);
            data.dataV = data.dataV.getSubset(1:data.nSamples);
            data.dataS = data.dataS.getSubset(1:data.nSamples);            
        end

        % no xsn bvh output
        data.dataX = [];

        % calib data
        
        qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
        data.calibS2B = mocapdb.XsensBody.loadCalib(data.calibFnameSensor2Body);
        data.calibS2B = data.calibS2B.adjustFrame(quatconj(qTCD2BM), [1 0 0 0], true);
        data.calibS2B.L_Foot.ori = quatmultiply(axang2quat([0 1 0 -pi/2]), ...
                            data.calibS2B.L_Foot.ori);
        data.calibS2B.R_Foot.ori = quatmultiply(axang2quat([0 1 0 -pi/2]), ...
                            data.calibS2B.R_Foot.ori);
        data.calibYawFix = struct();
        data.calibW2V = mocapdb.XsensBody.loadCalib(data.calibFnameSensorW2V);
        data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                           'w__x', zeros(1, 3));
                       
        % step data
        if exist(data.fnameRevStepDetect, 'file')
            data.revStepDetect = readtable(data.fnameRevStepDetect);
        else
            data.revStepDetect = table();
        end
        
        data.idx = [1, data.nSamples];
        
        if exist(dataDir, 'dir')
            if ~exist(sprintf('%s/mat', dataDir), 'dir')
                mkdir(sprintf('%s/mat', dataDir));
            end
            save(dataPath, 'data');
        end
    end
end

