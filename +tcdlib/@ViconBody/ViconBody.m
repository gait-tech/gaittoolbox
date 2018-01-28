% ======================================================================
%> @file TCDBody.m
%> @brief Body class for the TCD dataset
% ======================================================================

classdef ViconBody < mcViz.Body
    properties
        %> body is loaded from this source file name
        srcFileName
    end
    
    methods
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of TCDBody class.
        % ======================================================================
        function obj = ViconBody(varargin)
            obj@mcViz.Body(varargin{:});
        end
    end
end