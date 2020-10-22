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
        fs = 100
        % full trial start index
        ftStartIndex = 1
        % full trial end index
        ftEndIndex = inf
        
        % Head
        Head
        % Sternum
        Sternum
        % Pelvis
        Pelvis
        % Left upper arm
        L_UpArm
        % Right upper arm
        R_UpArm
        % Left low arm
        L_LowArm
        % Right low arm
        R_LowArm
        % Left thigh
        L_UpLeg
        % Right thigh
        R_UpLeg
        % Left shanks (near ankles)
        L_LowLeg
        % Left shanks (middle)
        L_LowLeg2
        % Right shanks (near ankles)
        R_LowLeg
        % Right shanks (middle)
        R_LowLeg2
        % Left foot
        L_Foot
        % Right foot
        R_Foot
    end
    
    properties (Constant)
		% segment property list
        segList = {'Head', 'Sternum', 'Pelvis', 'L_UpArm', 'R_UpArm', ...
            'L_LowArm', 'R_LowArm', 'L_UpLeg', 'R_UpLeg', ...
            'L_LowLeg', 'L_LowLeg2', 'R_LowLeg', 'R_LowLeg2', 'L_Foot', 'R_Foot'};
    end
    
    methods
        function obj = XsensBody(varargin)
            % Class constructor
            %
            % :param varargin: param1 (string), val1, param2 (string), val2
            %
            % :return: instance of BVHBody class.
            %
            % .. Author: - Luke Sy (UNSW GSBME)
        
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        function out = copyinfo(obj)
            out = mocapdb.XsensBody('srcFileName', obj.srcFileName, ...
                    'frame', obj.frame, 'nSamples', obj.nSamples, ...
                    'fs', obj.fs);
        end
        
        out = calcCalibSB(obj, refBody, sIdx);
        out = calcCalibSBFromMean(obj, refBody);
        out = calcCalibAnkleSensorW2PelvisWFromAcc(obj, idx);
        out = calcCalibAnkleSensorW2PelvisWFromROM(obj, calibS2B, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromGyroSkewness(obj, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromVicon(obj, dataV);
        out = exportRawMeasurementAsStruct(obj, seg, segAlias);
        out = getSubset(obj, idx);
        out = getSegSubset(obj, segList);
        out = getMean(obj);
        out = getNonemptySeg(obj);
        out = calcSegMeanRot(obj, seg, roty, idx);
        out = toViconFrame(obj, qR);
        out = adjustFrame(obj, qR1, qR2, orionly);
        out = calcGfrAcc(obj);
        out = conj(obj);
        out = removeAccBias(obj, bias);
        [out, bias] = normalizeAcc(obj, acc_norm);
        initializetoIdentity(obj);
        saveCalibCSV(obj, fname, mode);
        exportCSVs(obj, fname, info);
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
        [obj, idx] = loadCSVs(fname)
    end
end