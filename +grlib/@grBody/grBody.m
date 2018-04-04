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
    
    properties (Constant)
        posList = {'MIDPEL', 'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', ...
                   'RTIO'};
        oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK'};
    end
    
    methods (Hidden)
        function out = lim(obj, idx)
            low = inf; high = -inf;
            
            for i=1:length(obj.posList)
                data = obj.(obj.posList{i});
                low = min([low; data(:,idx)]);
                high = max([high;, data(:,idx)]);
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
        
        out = diff(obj1, obj2, seq);
        out = diffRMSE(obj1, obj2, seq);
    end
end