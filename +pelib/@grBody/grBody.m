% ======================================================================
%> @file Body.m
%> @brief Body class used to animate body and obtain gait parameters
% ======================================================================

classdef grBody < matlab.mixin.Copyable
    properties
        %> name of body
        name
        %> number of samples
        nSamples
        %> sampling frequency
        fs = 60
        %> frame: vicon / world / MIDPEL
        frame
        
        % Plot specifications
        posUnit = 'mm';
        oriUnit = 'deg';
        xyzColor = {'r', 'g', 'b'};
        axisScale = 0.25;
        rplColor = {'r', 'g', 'b'};
        lnSymbol = '-';
        ptSymbol = '.';
        
        %> SACR position (n x 3)
        MIDPEL
        LFEP
        LFEO
        LTIO
        LTOE
        RFEP
        RFEO
        RTIO
        RTOE
        
        %> pelvis orientation (n x 4)
        qRPV
        %> right femur orientation (n x 4)
        qRTH
        %> left femur orientation (n x 4)
        qLTH
        %> right tibia orientation (n x 4)
        qRSK
        %> left tibia orientation (n x 4)
        qLSK
        %> right foot orientation (n x 4)
        qRFT
        %> left foot orientation (n x 4)
        qLFT
    end
    
    properties (Constant)
        posList = {'MIDPEL', 'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', ...
                   'RTIO'};
        oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK'};
    end
    
    methods (Static)
        theta = calcJointAngles(prox, dist);
    end
    
    methods (Hidden)
        % ======================================================================
        %> @brief Return the minimum and maximum value of the corresponding coordinate
        %>
        %> @param idx 1 = x, 2 = y, 3 = z
        %>
        %> @return [low high]
        % ======================================================================
        function out = lim(obj, idx)
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
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of Body class.
        % ======================================================================
        function obj = grBody(varargin)
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
            if nargin <= 1, idx = 1:size(obj.qRPV, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRPV(idx, :), obj.qLTH(idx, :)) .* [-1 -1 -1];
        end
        
        function theta = calcJointAnglesRHip(obj, idx)
            if nargin <= 1, idx = 1:size(obj.qRPV, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRPV(idx, :), obj.qRTH(idx, :)) .* [1 -1 1];
        end

        function theta = calcJointAnglesLKnee(obj, idx)
            if nargin <= 1, idx = 1:size(obj.qLTH, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qLTH(idx, :), obj.qLSK(idx, :)) .* [-1 1 -1];
        end
        
        function theta = calcJointAnglesRKnee(obj, idx)
            if nargin <= 1, idx = 1:size(obj.qRTH, 1); end
            theta = pelib.grBody.calcJointAngles(obj.qRTH(idx, :), obj.qRSK(idx, :));
        end
        
        function d = getPelvisLength(obj, idx)
            if nargin <= 1, idx = 1; end
            d = norm(obj.LFEP(idx, :) - obj.RFEP(idx, :));
        end
        
        function d = getLShankLength(obj, idx)
            if nargin <= 1, idx = 1; end
            d = norm(obj.LFEO(idx, :) - obj.LTIO(idx, :));
        end
        
        function d = getRShankLength(obj, idx)
            if nargin <= 1, idx = 1; end
            d = norm(obj.RFEO(idx, :) - obj.RTIO(idx, :));
        end
        
        out = calcJointVel(obj, pts);
        out = calcJointAcc(obj, pts);
        out = calcDOri(obj, ref);
        
        out = diff(obj1, obj2, seq);
        out = diffRMSE(obj1, obj2, seq);
        out = toWorldFrame(obj, pos, ori);
        out = changePosUnit(obj, newUnit, update);
        out = getSubset(obj, idx);
        out = exportc3d(obj, fname, sensors, refBody, lsteps, rsteps, ...
                        extraMarkers, oriMode);
    end
end