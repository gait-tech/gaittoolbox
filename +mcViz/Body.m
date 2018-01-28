% ======================================================================
%> @file Body.m
%> @brief Body class used to animate body and obtain gait parameters
% ======================================================================

classdef Body < handle
    properties
        %> name of body
        name
        nSamples
        
        % Plot specifications
        posUnit = 'mm';
        oriUnit = 'deg';
        xyzColor = {'r', 'g', 'b'};
        lnSymbol = '-';
        ptSymbol = '.';
        
        %> SACR position (n x 3)
        SACR
        LFEP
        LFEO
        LTIO
        LTOE
        RFEP
        RFEO
        RTIO
        RTOE
        
        %> pelvis orientation (n x 4 OR 3 x 3 x n)
        qPelvis
        %> right femur orientation (n x 4 OR 3 x 3 x n)
        qRFemur
        %> left femur orientation (n x 4 OR 3 x 3 x n)
        qLFemur
        %> right tibia orientation (n x 4 OR 3 x 3 x n)
        qRTibia
        %> left tibia orientation (n x 4 OR 3 x 3 x n)
        qLTibia
        %> right foot orientation (n x 4 OR 3 x 3 x n)
        qRFoot
        %> left foot orientation (n x 4 OR 3 x 3 x n)
        qLFoot
    end
    
    properties (Constant)
        ptList = {'SACR', 'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', 'RTIO'};
    end
    
    methods (Hidden)
        function out = lim(obj, idx)
            low = inf; high = -inf;
            
            for i=1:length(obj.ptList)
                data = obj.(obj.ptList{i});
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
        function obj = Body(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
            if obj.SACR
                obj.nSamples = length(obj.SACR(:,1));
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
    end
    
    methods (Static)
        function plotPosComparison(varargin)
        end
        function plotPosDiff()
        end
        function plotOriComparison(varargin)
        end
        function plotOriDiff()
        end
    end
end