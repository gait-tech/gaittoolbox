classdef BVHBody < matlab.mixin.Copyable
	% a class for BVH body
    %
    % :param srcFileName: body is loaded from this source file name
	% :param frame: data are in this frame of reference (Vicon or IMU)
	
    properties
        srcFileName % body is loaded from this source file name
        frame % data are in this frame of reference (Vicon or IMU)
        posUnit = 'mm' % position unit
        nSamples % number of samples
        fs % sampling frequency
        
        Hips % similar to MIDPEL
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
        RightUpLeg % similar to RFEP
        RightLeg % similar to RFEO
        RightFoot % similar to RTIO
        RightToe % similar to RTOE
        LeftUpLeg % similar to LFEP
        LeftLeg % similar to LFEO
        LeftFoot % similar to LTIO
        LeftToe % similar to LTOE
        
        qHips % hip orientation (n x 4)
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
    
    properties (Hidden)
		% position property list
        posList = {'Hips', 'Spine', 'Spine1', 'Spine2', 'Spine3', ...
            'Neck', 'Head', ...
            'RightShoulder', 'RightArm', 'RightForeArm', 'RightHand', ...
            'LeftShoulder', 'LeftArm', 'LeftForeArm', 'LeftHand', ...
            'RightUpLeg', 'RightLeg', 'RightFoot', 'RightToe', ...
            'LeftUpLeg', 'LeftLeg', 'LeftFoot', 'LeftToe' };
		% orientation property list
        oriList = {'qHips', 'qSpine', 'qSpine1', 'qSpine2', 'qSpine3', ...
            'qNeck', 'qHead', ...
            'qRightShoulder', 'qRightArm', 'qRightForeArm', 'qRightHand', ...
            'qLeftShoulder', 'qLeftArm', 'qLeftForeArm', 'qLeftHand', ...
            'qRightUpLeg', 'qRightLeg', 'qRightFoot', 'qRightToe', ...
            'qLeftUpLeg', 'qLeftLeg', 'qLeftFoot', 'qLeftToe' };
    end
    
    methods
        % Class constructor
        %
        % :param varargin: param1 (string), val1, param2 (string), val2, ...
        %
        % :return: instance of BVHBody class.
		%
        % .. Author: - Luke Sy (UNSW GSBME)
		
        function obj = BVHBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
        
        out = togrBody(obj, idx, varargin)
        out = toViconBody(obj, idx, varargin)
    end
    
    methods (Hidden, Static)
        [skeleton,time] = loadbvh(fname)
    end
    
    methods (Static)
        obj = loadXsensBVHFile(fname, unit);
        obj = loadBVHFile(fname, unit);
        obj = loadOriPosFile(fname_ori, fname_pos, unit);
    end
    
    methods
        obj = toWorldFrame(obj, qR)
        out = getSubset(obj, idx)
        out = getStartIndex(obj)
        out = getEndIndex(obj)
        out = changePosUnit(obj, newUnit, update)
        out = adjustFootFrame2BMC(obj)
    end
end