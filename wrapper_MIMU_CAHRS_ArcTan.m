function [ varargout ] = wrapper_MIMU_CAHRS_ArcTan( varargin )
%WRAPPER_MIMU_CAHRS_ARCTAN Wrapper function designed to process a dataset of
%accelerometer, gyroscope and magnetometer through CAHRS
%   [ quatsMIMU ] = ...
%   wrapper_MIMU_CAHRS_ArcTan( muAcc,muMag,samplingRate, AccGyrMag )
%
%   Inputs:: (must be 'string', value pairs, e.g., 'fs',100 .... i.e., sampling rate of 100 Hz)
%       'qInit'          - initial orientation quaternion
%       'bDynamicMu'     - boolean vector length of the recording, i.e., length(acc(:,1))
%       'static_mu_acc'  - \mu_{a} which will be used when 'bDynamicMu' is false
%       'static_mu_mag'  - \mu_{m} which will be used when 'bDynamicMu' is false
%       'dynamic_mu_acc' - \mu_{a} which will be used when 'bDynamicMu' is true
%       'dynamic_mu_mag' - \mu_{a} which will be used when 'bDynamicMu' is true
%       'fs'             - sampling rate
%       'Acc'            - accelerometer signal size = [N x 3]
%       'Gyr'            - gyroscope     signal size = [N x 3]
%       'Mag'            - magnetometer  signal size = [N x 3]
%
%   Outputs::
%       quatsMIMU - quaternions from a 9DOF MIMU

qInit = [];
boolDynamicMu = [];
for i = 1:2:nargin
    if strcmp(varargin{i}, 'qInit'),
        qInit = varargin{i+1};
    elseif strcmp(varargin{i},'static_mu_acc'),  
        static_mu_acc = varargin{i+1};      
    elseif strcmp(varargin{i},'static_mu_mag'), 
        static_mu_mag = varargin{i+1};
    elseif strcmp(varargin{i},'dynamic_mu_acc'),  
        dynamic_mu_acc = varargin{i+1};      
    elseif strcmp(varargin{i},'dynamic_mu_mag'), 
        dynamic_mu_mag = varargin{i+1};        
    elseif strcmp(varargin{i},'fs'), 
        fs = varargin{i+1};   
    elseif strcmp(varargin{i},'Acc')
        ACC = varargin{i+1};
        [NUM_ACC,COLS_A]=size(ACC);
        if COLS_A ~= 3
            error('Accelerometer data must be size [N x 3]');
        end        
    elseif strcmp(varargin{i},'Gyr')
        GYRO = varargin{i+1};
        [NUM_GYR,COLS_G]=size(GYRO);
        if COLS_G ~= 3
            error('Gyroscope data must be size [N x 3]');
        end
    elseif strcmp(varargin{i},'Mag')
        MAGNO = varargin{i+1};        
        [NUM_MAG,COLS_M]=size(MAGNO);
        if COLS_M ~= 3
            error('Magnetometer data must be size [N x 3]');
        end
    elseif strcmp(varargin{i},'bDynamicMu');
        boolDynamicMu = varargin{i+1};
    else error('Invalid argument');
    end    
end

if (NUM_ACC ~= NUM_GYR) && (NUM_ACC ~= NUM_MAG) && (NUM_GYR ~= NUM_MAG)
    error('Number of samples in IMU data are unequal');
end

[FRAMES,~] = size(ACC);
quatsMIMU  = zeros(FRAMES,4);

%% Dynamic Tuning
if ~isempty(boolDynamicMu)
    cahrs = CAHRS('fs',fs,'muAcc',0.5,'muMag',0.5,'q0',qInit);
    for frameNum = 1:FRAMES
        % if true, dyamic - use the specified muAcc and muMag
        if boolDynamicMu(frameNum)
            cahrs.muAcc = dynamic_mu_acc;      cahrs.muMag = dynamic_mu_mag;
        % if false, use the static_mu_acc and static_mu_mag
        else
            cahrs.muAcc = static_mu_acc;      cahrs.muMag = static_mu_mag;
        end
        cahrs.UpdateNoCallsArcTan(GYRO(frameNum,:), ACC(frameNum,:),MAGNO(frameNum,:));
        quatsMIMU(frameNum, :) = cahrs.qGlobal;
    end
else
    error('Boolean Vector must be specified');
end
% variable output assignment
if nargout >= 1
    varargout{1} = quatsMIMU;
end
% add other output arguments here
end




