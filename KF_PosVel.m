classdef KF_PosVel < handle
    %KF_POSVEL class definition for kalman filter to estimate position and
    %velocity incorporating traditional kalman filtering and a third state
    %in which state constraints are enforced
    %   Detailed explanation goes here
    
    % Class members here are not changeable unless get/set methods are
    % declared
    properties(Access = private)    
        A = eye(2);     % state transition matrix - diagonal 6x6 matrix
        B = zeros(2,1); % input matrix
    end

    % Class members here can be changed at any time
    properties(Access = public)
        % sampling rate
        fs = NaN;
        % inter-arrival time
        dt = NaN;          
        % state vector, position & velocity in R3
        xHat = zeros(2,1); 
        % variance in the accelerometer measurement which is the input to
        % the KF
        sigma_acc = NaN;   
        Q = NaN;
        P = NaN;
        H = zeros(1,2);    % measurement matrix
        y = NaN;
        Rvel = NaN;
        S = NaN;
        K = NaN;
    end
    
    methods (Access = public)
        function obj = KF_PosVel(varargin)
            for i = 1:2:nargin
                if     strcmp('fs',varargin{i}),
                    obj.fs = varargin{i+1};
                    obj.dt = 1/obj.fs;
                    obj.A(1,2) = obj.dt;
                    obj.B(1,1) = 0.5 * obj.dt^2;
                    obj.B(2,1) = obj.dt;
                elseif strcmp('x0',varargin{i}),
                    obj.xHat = varargin{i+1};                                        
                elseif strcmp('sigma_acc',varargin{i}),
                    obj.sigma_acc = varargin{i+1};
                else
                    error('Not Supported');
                end
            end
            % calculate the process noise covariance matrix
            obj.Q = obj.B * obj.B'* obj.sigma_acc^2;
            obj.P = obj.Q;
            obj.H(1,2) = 1; 
            if isnan(obj.sigma_acc), error('sigma_acc undefined'); end
        end
        % Assuming a fixed sampling rate
        function obj = PredictStep(obj, acc)
            [Racc,Cacc] = size(acc);
            if ~isequal(Racc * Cacc,1), error('Input acceleration must be 3 x 1 vector');
            else if Cacc == 1, acc=acc'; end
            end
            obj.xHat = obj.A * obj.xHat + obj.B * acc;
            obj.P    = obj.A * obj.P * obj.A' + obj.Q;
        end
        % Assuming a fixed sampling rate, perform a zero velcity update,
        % asumming a noisy zero value for the measured velocity
        function obj = UpdateStep(obj, sigma_vel)
            obj.y = 0 - obj.xHat(2); % innnovation/ measurement residual
            obj.Rvel = sigma_vel^2;
            obj.S = obj.H * obj.P * obj.H' + obj.Rvel; % scalar
            obj.K = obj.P * obj.H' * obj.S^(-1);
            obj.xHat = obj.xHat + obj.K * obj.y;
            obj.P    = (eye(2)-obj.K*obj.H)*obj.P;
        end        
    end
    
end

