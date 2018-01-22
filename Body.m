% ======================================================================
%> @file Body.m
%> @brief Body class used to animate body and obtain gait parameters
% ======================================================================

classdef Body < handle
    properties
        %> name of body
        name
        nSamples
        
        %> SACR position (n x 3)
        SACR
        LFEP
        LFEO
        LTIO
        RFEP
        RFEO
        RTIO
        
        %> pelvis orientation (n x 4 OR n x 3 x 3)
        qPelvis
        %> right femur orientation (n x 4 OR n x 3 x 3)
        qRFemur
        %> left femur orientation (n x 4 OR n x 3 x 3)
        qLFemur
        %> right tibia orientation (n x 4 OR n x 3 x 3)
        qRTibia
        %> left tibia orientation (n x 4 OR n x 3 x 3)
        qLTibia
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
            
            obj.nSamples = length(obj.SACR(:,1));
        end
        
        function calcJointAngles()
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