function data = loadNeuRaTrialData(dataDir, subjName, actName, ns, usebuffer)
    % Load data for selected movement trial
    % 
    % :param dataDir: data directory path
    % :param subjName: subject name of trial movement
    % :param actName: action name of trial movement
    % :param ns: yaw offset calibration mode. 
    %            NS1 = use yaw offset from ROM
    %            NS2 = use yaw offset from Vicon
    % :param usebuffer: [optional] use buffer .mat if available
    % :param options: [optional] id to body mapping. By default uses the 
    %                 map specified in https://gait-tech.github.io/gaittoolbox/data.html.
    % 
    % :return: struct data with trial data inside
    %
    % .. Author: - Luke Wicent Sy (GSBME, Modified 2020 May 8)
    if nargin <= 3
        ns = '';
    end
    if nargin <= 4
        usebuffer = true;
    end
    
    if strcmp(ns, "")
        name = sprintf("%s-%s", subjName, actName);
    else
        name = sprintf("%s-%s-%s", ns, subjName, actName);
    end
    
    dataPath = sprintf('%s/mat/%s-%s-%s-buffer.mat', dataDir, ...
                        extractBetween(ns,1,3), subjName, actName);
    if usebuffer && exist(dataPath, 'file')
        load(dataPath, 'data');
    else
        data = struct('name', name, ...
            'fnameV', sprintf('%s/vicon/%s-%s', dataDir, subjName, actName), ...
            'fnameX', sprintf('%s/xsens/%s-%s.bvh', dataDir, subjName, actName), ...
            'fnameS', sprintf('%s/imu/%s-%s', dataDir, subjName, actName), ...
            'calibFnameSensorYawFixWorldFrame', ...
            sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrame.txt', ...
                    dataDir, subjName, actName), ...
            'calibFnameSensorW2V', ...
            sprintf('%s/calib/%s-Calib-SensorW2V.txt', dataDir, subjName), ...
            'fnameRevStepDetect', ...
            sprintf('%s/rawstep-detect/%s-%s-revStepDetect.csv', ...
                    dataDir, subjName, actName));

        data.dataV = mocapdb.ViconBody.loadCSV(data.fnameV); 
        % load imu measurements (acc, gyr, mag are in sensor frame, 
        % ori = orientation of sensor in world frame)
        [data.dataS, data.idx] = mocapdb.XsensBody.loadCSVs(data.fnameS); 

        data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
        % transformation to express dataX in biomechanical convention
        % x = forward, z = upward, y = left
        qXsens2BMC = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        data.dataX = dataX.getSubset(data.idx(1):data.idx(2)) ...
                            .toWorldFrame(qXsens2BMC) ...
                            .changePosUnit('m');

        data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
        data.calibW2V = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorW2V);
        data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                           'w__x', zeros(1, 3));
        data.revStepDetect = readtable(data.fnameRevStepDetect);
    end
end

