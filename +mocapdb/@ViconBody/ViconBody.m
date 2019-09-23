% ======================================================================
%> @file ViconBody.m
%> @brief Body class for Vicon csv export
% ======================================================================

classdef ViconBody < matlab.mixin.Copyable
    properties
        %> body is loaded from this source file name
        srcFileName
        %> data are in this frame of reference (Vicon or IMU)
        frame
        %> sampling rate
        fs
        posUnit = 'mm'
        %> number of samples
        nSamples
        
        %> Position
        PELV
        LFEP
        LFEO
        LTIO
        LTOE
        RFEP
        RFEO
        RTIO
        RTOE
        
        %> Orientation
        %> pelvis orientation (n x 4 OR 3 x 3 x n)
        qRPV
        %> right femur orientation (n x 4 OR 3 x 3 x n)
        qRTH
        %> left femur orientation (n x 4 OR 3 x 3 x n)
        qLTH
        %> right tibia orientation (n x 4 OR 3 x 3 x n)
        qRSK
        %> left tibia orientation (n x 4 OR 3 x 3 x n)
        qLSK
        %> right foot orientation (n x 4 OR 3 x 3 x n)
        qRFT
        %> left foot orientation (n x 4 OR 3 x 3 x n)
        qLFT
    end
    
    properties (Hidden)
        posList = {'PELV', 'LFEP', 'LFEO', 'LTIO', 'LTOE', ...
            'RFEP', 'RFEO', 'RTIO', 'RTOE'};
        oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK'};
    end
    methods
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of ViconBody class.
        % ======================================================================
        function obj = ViconBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        out = togrBody(obj, idx, args);
        startIdx = getStartIndex(obj);
        endIdx = getEndIndex(obj);
        out = getSubset(obj, idx);
        out = changePosUnit(obj, newUnit, update);
    end
    
    methods (Static)
        obj = loadViconMat(fname)
    end
end