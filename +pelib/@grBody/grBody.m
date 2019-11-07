classdef grBody < matlab.mixin.Copyable
	% Body class used to animate body and obtain gait parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    properties
        % name of body
        name
        % number of samples
        nSamples
        % sampling frequency
        fs = 60
        % frame: vicon / world / MIDPEL
        frame
        
        % Plot specifications
        posUnit = 'mm';
        oriUnit = 'deg';
        xyzColor = {'r', 'g', 'b'};
        axisScale = 0.25;
        rplColor = {'r', 'g', 'b'};
        lnSymbol = '-';
        ptSymbol = '.';
        
        % SACR position (n x 3)
        MIDPEL
        % Left hip joint center
        LFEP
        % Left knee joint center
        LFEO
        % Left ankle joint center
        LTIO
        % Left toe
        LTOE
        % Right hip joint center
        RFEP
        % Right knee joint center
        RFEO
        % Right ankle joint center
        RTIO
        % Right toe
        RTOE
        
        % pelvis orientation (n x 4), z = upward, x = forward, y = towards left hip joint center
        qRPV
        % right femur orientation (n x 4), z = along thigh, x = forward, y towards left
        qRTH
        % left femur orientation (n x 4), z = along thigh, x = forward, y towards left
        qLTH
        % right tibia orientation (n x 4), z = along tibia, x = forward, y towards left
        qRSK
        % left tibia orientation (n x 4), z = along tibia, x = forward, y towards left
        qLSK
        % right foot orientation (n x 4) following Vicon convention, 
        % z = toe to heel vector projected to the floor, x = upward, y = towards left
        qRFT
        % left foot orientation (n x 4) following Vicon convention
        % z = toe to heel vector projected to the floor, x = upward, y = towards left
        qLFT
    end
    
    properties (Constant)
		% position property list
        posList = {'MIDPEL', 'LFEP', 'LFEO', 'LTIO', 'LTOE', ...
                   'RFEP', 'RFEO', 'RTIO', 'RTOE'};
	    % orientation property list
        oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK', 'qRFT', 'qLFT'};
    end
    
    methods (Static)
        theta = calcJointAngles(prox, dist);
        dist = calcDistRotm(prox, angles, seq);
        prox = calcProxRotm(dist, angles, seq);
        out = generateBodyFromJointAngles(posMP, qOriMP, ...
                        anglesLT, anglesRT, angleLK, angleRK, ...
                        dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, seq);
        [alphaLK, alphaRK] = calcKneeAnglesFromMPLARADist( ...
                        PELV_CS, LTIB_CS, RTIB_CS, ...
                        dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, ...
                        dMPLADist, dMPRADist);
    end
    
    methods (Hidden)	
        function out = lim(obj, idx)
			% Return the minimum and maximum value of the corresponding coordinate
			%
			% :param idx: 1 = x, 2 = y, 3 = z
			%
			% :return: [low high]
			%
			% .. Author: - Luke Sy (UNSW GSBME)
            low = inf; high = -inf;
            
            for i=1:length(obj.posList)
                data = obj.(obj.posList{i});
                low = min([low; data(:,idx)]);
                high = max([high; data(:,idx)]);
            end
            
            out = [low high];
        end
    end
    
    methods
        function obj = grBody(varargin)
			% Class constructor
			%
			% :param varargin: param1 (string), val1, param2 (string), val2, ...
			%
			% :return: instance of Body class.
			%
			% .. Author: - Luke Sy (UNSW GSBME)
		
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
            if obj.MIDPEL
                obj.nSamples = length(obj.MIDPEL(:,1));
            end
        end
        
        function out = xlim(obj)
            out = lim(obj, 1);
        end
        
        function out = ylim(obj)
            out = lim(obj, 2);
        end
        
        function out = zlim(obj)
            out = lim(obj, 3);
        end
        
        function [ptsX, ptsY, ptsZ] = groundCoordinates(obj)
%             xVector = obj.RTIO(1,:) - obj.LTIO(1,:);
%             xVector = xVector / norm(xVector);
            gndOrigin = (obj.LTIO(1,:)+obj.RTIO(1,:))/2;
%             zVector = obj.MIDPEL(1,:) - gndOrigin;
%             yVector = cross(zVector, xVector);
            CS = quat2rotm(obj.qRPV(1,:));
            xVector = CS(:,1)';
            yVector = CS(:,2)';
            pts = [ gndOrigin + xVector + yVector;
                    gndOrigin - xVector + yVector;
                    gndOrigin - xVector - yVector;
                    gndOrigin + xVector - yVector];
            
            ptsX = pts(:,1); ptsY = pts(:,2); ptsZ = pts(:,3);
        end

        function theta = calcJointAnglesLHip(obj, idx)
			% Calculate left hip joint angles (seq: YXZ)
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: all)
			%
			% :return: theta - joint angles (n x 3)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1:size(obj.qRPV, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRPV(idx, :), obj.qLTH(idx, :));
            theta = theta(:, [2 1 3]) .* [-1 -1 -1];
        end
        
        function theta = calcJointAnglesRHip(obj, idx)
			% Calculate right hip joint angles (seq: YXZ)
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: all)
			%
			% :return: theta - joint angles (n x 3)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1:size(obj.qRPV, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRPV(idx, :), obj.qRTH(idx, :));
            theta = theta(:, [2 1 3]) .* [1 -1 1];
        end

        function theta = calcJointAnglesLKnee(obj, idx)
			% Calculate left knee joint angles (seq: YXZ)
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: all)
			%
			% :return: theta - joint angles (n x 3)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1:size(obj.qLTH, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qLTH(idx, :), obj.qLSK(idx, :));
            theta = theta(:, [2 1 3]) .* [-1 1 -1];
        end
        
        function theta = calcJointAnglesRKnee(obj, idx)
			% Calculate right knee joint angles (seq: YXZ)
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: all)
			%
			% :return: theta - joint angles (n x 3)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1:size(obj.qRTH, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRTH(idx, :), obj.qRSK(idx, :));
            theta = theta(:, [2 1 3]);
        end
        
        function d = getPelvisLength(obj, idx)
			% Calculate pelvis length
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: 1)
			%
			% :return: d - body segment length (n x 1)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1; end
            d = norm(obj.LFEP(idx, :) - obj.RFEP(idx, :));
        end
        
        function d = getLShankLength(obj, idx)
			% Calculate left shank length
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: 1)
			%
			% :return: d - body segment length (n x 1)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1; end
            d = norm(obj.LFEO(idx, :) - obj.LTIO(idx, :));
        end
        
        function d = getRShankLength(obj, idx)
			% Calculate right shank length
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: 1)
			%
			% :return: d - body segment length (n x 1)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            if nargin <= 1, idx = 1; end
            d = norm(obj.RFEO(idx, :) - obj.RTIO(idx, :));
        end
        
        function d = calcLFemurLength(obj)
			% Calculate left femur length
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: 1)
			%
			% :return: d - body segment length (n x 1)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            d = vecnorm(obj.LFEP - obj.LFEO, 2, 2);
        end
        
        function d = calcRFemurLength(obj)
			% Calculate right femur length
			%
			% :param obj: this object class
			% :param idx: [OPTIONAL] index to be calculated (default: 1)
			%
			% :return: d - body segment length (n x 1)
			%
			% .. Author: - Luke Sy (UNSW GSBME)
			
            d = vecnorm(obj.RFEP - obj.RFEO, 2, 2);
        end
        
        out = calcJointVel(obj, pts);
        out = calcSegAngVel(obj, pts, frame);
        out = calcJointAcc(obj, pts);
        out = calcDOri(obj, ref);
        [out1, out2] = calcDOrinobias(obj, ref);
        out = calcDPos(obj, ref);
        out = calcMPLARAdist(obj);
        
        out = diff(obj1, obj2, seq);
        out = diffRMSE(obj1, obj2, seq);
        out = diffRMSEandMean(obj1, obj2);
        out = toWorldFrame(obj, pos, ori);
        out = changePosUnit(obj, newUnit, update);
		out = changeRefFrame(obj, ref)
        out = getSubset(obj, idx);
        out = exportc3d(obj, fname, sensors, refBody, lsteps, rsteps, ...
                        extraMarkers, oriMode, spevents);
    end
end