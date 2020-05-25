classdef XsensBody < matlab.mixin.Copyable
	% Xsens Body class for the TCD dataset
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    properties
        % body is loaded from this source file name
        srcFileName
        % data are in this frame of reference (Vicon or World or Calib)
        frame
        % number of samples
        nSamples
        % sampling frequency
        fs = 100;
        
        Head
        Sternum
        Pelvis
        L_UpArm
        R_UpArm
        L_LowArm
        R_LowArm
        L_UpLeg
        R_UpLeg
        L_LowLeg
<<<<<<< HEAD
        R_LowLeg
=======
        L_LowLeg2
        R_LowLeg
        R_LowLeg2
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
        L_Foot
        R_Foot
    end
    
    properties (Constant)
		% segment property list
        segList = {'Head', 'Sternum', 'Pelvis', 'L_UpArm', 'R_UpArm', ...
            'L_LowArm', 'R_LowArm', 'L_UpLeg', 'R_UpLeg', ...
<<<<<<< HEAD
            'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
=======
            'L_LowLeg', 'L_LowLeg2', 'R_LowLeg', 'R_LowLeg2', 'L_Foot', 'R_Foot'};
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
    end
    
    methods
        % Class constructor
        %
        % :param varargin: param1 (string), val1, param2 (string), val2, ...
        %
        % :return: instance of BVHBody class.
        %
		% .. Author: - Luke Sy (UNSW GSBME)
        function obj = XsensBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
<<<<<<< HEAD
        out = calcCalibSB(obj, refBody, sIdx);
=======
        function out = copyinfo(obj)
            out = mocapdb.XsensBody('srcFileName', obj.srcFileName, ...
                    'frame', obj.frame, 'nSamples', obj.nSamples, ...
                    'fs', obj.fs);
        end
        
        out = calcCalibSB(obj, refBody, sIdx);
        out = calcCalibSBFromMean(obj, refBody);
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
        out = calcCalibAnkleSensorW2PelvisWFromAcc(obj, idx);
        out = calcCalibAnkleSensorW2PelvisWFromROM(obj, calibS2B, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromGyroSkewness(obj, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromVicon(obj, dataV);
        out = exportRawMeasurementAsStruct(obj, seg, segAlias);
        out = getSubset(obj, idx);
<<<<<<< HEAD
        out = toViconFrame(obj, qR);
        initializetoIdentity(obj);
        saveCalibCSV(obj, fname);
=======
        out = getSegSubset(obj, segList);
        out = toViconFrame(obj, qR);
        out = adjustFrame(obj, qR1, qR2, orionly);
        out = calcGfrAcc(obj);
        out = conj(obj);
        initializetoIdentity(obj);
        saveCalibCSV(obj, fname);
        exportCSVs(obj, fname, info)
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
    end
    
    methods (Hidden, Static)
        obj = load_mvnx(fname)
    end
    
    methods (Static)
        obj = loadSensorFile(fname)
        obj = loadCalib(fname)
        obj = loadMTExport(name, options)
        obj = loadMVNX(fname, options)
        obj = loadCalibSensorW2V(viconFName, xsensFName, options, idx)
        obj = loadCalibCSV(obj, fname)
<<<<<<< HEAD
=======
        [obj, idx] = loadCSVs(fname)
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
    end
end