% ======================================================================
%> @file TCDBody.m
%> @brief Body class for the TCD dataset
% ======================================================================

classdef BVHBody < handle
    properties
        %> body is loaded from this source file name
        srcFileName
        %> data are in this frame of reference (Vicon or IMU)
        frame
        posUnit = 'mm'
        %> number of samples
        nSamples
        %> sampling frequency
        fs
        
        %> similar to MIDPEL
        Hips
        Spine
        Spine1
        Spine2
        Spine3
        Neck
        Head
        RightShoulder
        RightArm
        RightForeArm
        RightHand
        LeftShoulder
        LeftArm
        LeftForeArm
        LeftHand
        %> similar to RFEP
        RightUpLeg
        %> similar to RFEO
        RightLeg
        %> similar to RTIO
        RightFoot
        RightToe
        %> similar to LFEP
        LeftUpLeg
        %> similar to LFEO
        LeftLeg
        %> similar to LTIO
        LeftFoot
        LeftToe
        
        %> hip orientation (n x 4)
        qHips
        qSpine
        qSpine1
        qSpine2
        qSpine3
        qNeck
        qHead
        qRightShoulder
        qRightArm
        qRightForeArm
        qRightHand
        qLeftShoulder
        qLeftArm
        qLeftForeArm
        qLeftHand
        qRightUpLeg
        qRightLeg
        qRightFoot
        qRightToe
        qLeftUpLeg
        qLeftLeg
        qLeftFoot
        qLeftToe
    end
    
    methods
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of BVHBody class.
        % ======================================================================
        function obj = BVHBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        out = togrBody(obj, varargin)
    end
    
    methods (Hidden, Static)
        [skeleton,time] = loadbvh(fname)
    end
    
    methods (Static)
        obj = loadXsensBVHFile(fname, unit)
        obj = loadBVHFile(fname, unit)
        obj = loadOriPosFile(fname_ori, fname_pos, unit)
    end
    
    methods
        obj = toWorldFrame(obj, qR)
    end
end