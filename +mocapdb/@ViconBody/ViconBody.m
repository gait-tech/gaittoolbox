classdef ViconBody < matlab.mixin.Copyable
	% Body class from Vicon csv export
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    properties
        srcFileName % body is loaded from this source file name
        frame % data are in this frame of reference (Vicon or IMU)
        fs % sampling rate
        posUnit = 'mm' % position unit
        nSamples % number of samples
        ftStartIndex = 1 % full trial start index
        ftEndIndex = inf % full trial end index
        
        %% Position:
        PELV % pelvis position
        LFEP % left hip position
        LFEO % left knee position
        LTIO % left ankle position
        LTOE % left toe position
        RFEP % right hip position
        RFEO % right knee position
        RTIO % right ankle position
        RTOE % right toe position
        
        %% Orientation:
        qRPV % pelvis orientation (n x 4 OR 3 x 3 x n)
        qRTH % right femur orientation (n x 4 OR 3 x 3 x n)
        qLTH % left femur orientation (n x 4 OR 3 x 3 x n)
        qRSK % right tibia orientation (n x 4 OR 3 x 3 x n)
        qLSK % left tibia orientation (n x 4 OR 3 x 3 x n)
        qRFT % right foot orientation (n x 4 OR 3 x 3 x n)
        qLFT % left foot orientation (n x 4 OR 3 x 3 x n)
    end
    
    properties (Hidden)
		% position property list
        posList = {'PELV', 'LFEP', 'LFEO', 'LTIO', 'LTOE', ...
            'RFEP', 'RFEO', 'RTIO', 'RTOE'};
		% orientation property list
        oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK', 'qRFT', 'qLFT'};
    end
    
    methods
        % Class constructor
        %
        % :param varargin: param1 (string), val1, param2 (string), val2, ...
        %
        % :return: instance of ViconBody class.
		%
        % .. Author: - Luke Sy (UNSW GSBME)
		
        function obj = ViconBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        out = togrBody(obj, idx, args);
        startIdx = getStartIndex(obj);
        endIdx = getEndIndex(obj, untilfirstnan);
        out = getSubset(obj, idx);
        out = changePosUnit(obj, newUnit, update);
        exportCSV(obj, fname, info);
    end
    
    methods (Static)
        obj = loadViconMat(fname, plugingait)
        [obj, idx] = loadCSV(fname)
    end
end