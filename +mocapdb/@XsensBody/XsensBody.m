% ======================================================================
%> @file XsensBody.m
%> @brief Xsens Body class for the TCD dataset
% ======================================================================

classdef XsensBody < matlab.mixin.Copyable
    properties
        %> body is loaded from this source file name
        srcFileName
        %> data are in this frame of reference (Vicon or World or Calib)
        frame
        %> number of samples
        nSamples
        %> sampling frequency
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
        R_LowLeg
        L_Foot
        R_Foot
    end
    
    properties (Constant)
        segList = {'Head', 'Sternum', 'Pelvis', 'L_UpArm', 'R_UpArm', ...
            'L_LowArm', 'R_LowArm', 'L_UpLeg', 'R_UpLeg', ...
            'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    end
    
    methods
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of BVHBody class.
        % ======================================================================
        function obj = XsensBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        out = calcCalibSB(obj, refBody, sIdx);
        out = calcCalibAnkleSensorW2PelvisWFromAcc(obj, idx);
        out = calcCalibAnkleSensorW2PelvisWFromROM(obj, calibS2B, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromGyroSkewness(obj, DEGRANGE);
        out = calcCalibAnkleSensorW2PelvisWFromVicon(obj, dataV);
        out = exportRawMeasurementAsStruct(obj, seg, segAlias);
        out = getSubset(obj, idx);
        out = toViconFrame(obj, qR);
        initializetoIdentity(obj);
        saveCalibCSV(obj, fname);
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
    end
end