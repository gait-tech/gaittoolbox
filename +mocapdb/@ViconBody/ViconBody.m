classdef ViconBody < matlab.mixin.Copyable
	% Body class from Vicon csv export
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    properties
        % body is loaded from this source file name
        srcFileName
        % data are in this frame of reference (Vicon or IMU)
        frame
        % sampling rate
        fs
		% position unit
        posUnit = 'mm'
        % number of samples
        nSamples
        % full trial start index
        ftStartIndex = 1
        % full trial end index
        ftEndIndex = inf
        
		% Positions:
        % pelvis position
        PELV
		% left hip position
        LFEP
		% left knee position
        LFEO
		% left ankle position
        LTIO
		% left toe position
        LTOE
		% right hip position
        RFEP
		% right knee position
        RFEO
		% right ankle position
        RTIO
		% right toe position
        RTOE
        
        % Orientation:
        % pelvis orientation (n x 4 OR 3 x 3 x n)
        qRPV
        % right femur orientation (n x 4 OR 3 x 3 x n)
        qRTH
        % left femur orientation (n x 4 OR 3 x 3 x n)
        qLTH
        % right tibia orientation (n x 4 OR 3 x 3 x n)
        qRSK
        % left tibia orientation (n x 4 OR 3 x 3 x n)
        qLSK
        % right foot orientation (n x 4 OR 3 x 3 x n)
        qRFT
        % left foot orientation (n x 4 OR 3 x 3 x n)
        qLFT
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