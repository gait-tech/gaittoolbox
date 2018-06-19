% ======================================================================
%> @file XsensBody.m
%> @brief Xsens Body class for the TCD dataset
% ======================================================================

classdef XsensBody < handle
    properties
        %> body is loaded from this source file name
        srcFileName
        %> data are in this frame of reference (Vicon or World or Calib)
        frame
        %> number of samples
        nSamples
        
        Head
        Sternum
        Pelvis
        L_UpArm
        R_UpArm
        L_LowArm
        R_LowArm
        L_UpLeg
        R_UpLeg
        L_LowLeg
        R_LowLeg
        L_Foot
        R_Foot
    end
    
    methods
        % ======================================================================
        %> @brief Class constructor
        %>
        %> @param varargin param1 (string), val1, param2 (string), val2, ...
        %>
        %> @return instance of BVHBody class.
        % ======================================================================
        function obj = XsensBody(varargin)
            for i = 1:2:nargin
               obj.(varargin{i}) = varargin{i+1};
            end
        end
    end
    
    methods (Hidden, Static)
    end
    
    methods (Static)
        obj = loadSensorFile(fname)
        obj = loadCalib(fname)
        obj = loadMTExport(name, options)
    end
end