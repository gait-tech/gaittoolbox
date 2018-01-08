classdef CAHRS <handle
%CAHRS is the AHRS algorithm for a 9DOF and or 6DOF MIMU/IMU 
% See 'Quaternion based Complementary Filter for Attitude Determination of a
% Smartphone'
% information from each of the sensors has been decoupled and fused
% according to the tuning parameters specified by the user. The
% algorithm assumes that all sensors inputs have been calibrated
%
% Author(s): Stephen Redmond ()
%            Michael Del Rosario (m.delrosario@unsw.edu.au)
%            Phils Ngo - developed/tested Update_IKF

properties (Access = private)    
    ctr       = 0; % idx corresponding to the number of samples processed
    sigma_gyr = NaN;
    sigma_acc = NaN;
    R_gamma_a_lower = NaN;
    
    new_avg_gam_a = 0;
    pre_avg_gam_a = 0;
    var_sum_gam_a = 0;
    B_VAR_GAM_A = [];    
    N_VAR_GAM_A = NaN;    

    new_avg_gam_m = 0;
    pre_avg_gam_m = 0;
    var_sum_gam_m = 0;    
    B_VAR_GAM_M = [];
    N_VAR_GAM_M = NaN;  
    
    gRef = 9.81;
    c_a = 0.1;
end


properties (Access = public)
    fs          = 100; % additional constants defined for speed
    dt          = 0.01;
    qAngVelDev  = [1 0 0 0];
    qGlobal     = [1 0 0 0];
    qGlobalPrev = [1 0 0 0];        
    muAcc = 0.005;
    muMag = 0.005;
    qGlobalGyr = NaN;
    AccGlobal = NaN;
    qRotTowardsUp = NaN;
    MagGlobal = NaN;
    qRotTowardsNorth = NaN;
    qGlobalPrevReal  = NaN;
    qGlobalRotatedUp = NaN;
    
    % these are for storing the angles gamma_a and gamma_m
    gammaA         = NaN;
    gammaM         = NaN;    
    cos_muA_gammaA = NaN;
    sin_muA_gammaA = NaN;
    cos_muM_gammaM = NaN;
    sin_muM_gammaM = NaN;    
    % store global acceleration vector
    gAcc = [];
    %% New variables for Adaptive gain via indirect kalman filter
    % --- Process Noise Covariance Matrix, in this case 1D --- %
    % Q = 2 * dt^{2} * \sigma_{gyr}^{2} i.e., the variance in the gyro when
    % stationary, user should specify the noise level of the gyro for
    % example, sigma = 0.015
    Q      = NaN; 
    R_gamma_a = NaN; % variance in the measurement noise  
    avg_gamma_a = NaN;
    P_pri  = NaN; % A Priori error state covariance matrix
    P_pos  = NaN; % A Posteriori error state covariance matrix, can be initialised to Q
    mu_a   = NaN; % adaptive kalman gain
    
    cos_K_gammaA = NaN;    
    sin_K_gammaA = NaN;    
    avg_gam_a = NaN;
    var_gam_a = NaN;

    Qmag = NaN; 
    R_gamma_m = NaN; % variance in the measurement noise 
    avg_gamma_m = NaN;
    Pmag_pri = NaN; % initialise a posteriori as Q        
    Pmag_pos = NaN; % initialise a priori as Q
    mu_m = NaN;

    cos_K_gammaM = NaN;    
    sin_K_gammaM = NaN; 
    avg_gam_m = NaN;
    var_gam_m = NaN;    

%     N_VAR_GAM_M = 100; 
%     B_VAR_GAM_M = nan(N_VAR_GAM_M,1); % weighted average for gamma_m
%     avg_gam_m = 1/N_VAR_GAM_M;        % weighted average for variance of gamma_m
%     var_gam_m = 1/(N_VAR_GAM_M-1);    
    ext_a = [0, 0, 0];    
end

methods (Access = public)
    function obj = CAHRS(varargin)
        for i = 1:2:nargin
            if strcmp(varargin{i},'fs'), 
                obj.fs = varargin{i+1};
                obj.dt=1/varargin{i+1};
            elseif strcmp(varargin{i},'muAcc'),
                obj.muAcc = varargin{i+1};                
            elseif strcmp(varargin{i},'muMag'), 
                obj.muMag = varargin{i+1};
            elseif strcmp(varargin{i},'q0'),
                obj.qGlobalPrev = varargin{i+1};
                obj.qGlobal = varargin{i+1};
%% --- CAHRS_Indirect Kalman Filter ---   
% add additional fields to support Indirect Kalman Filter
% constructor for CAHRS with indirect kalman filter for dynamic correction
% rate of mu_a and mu_m based on the variance of gamma_a and gamma_m                
            elseif strcmp(varargin{i},'N_VAR_GAM_A');
                % number of samples to calculate variance in gamma_a,
                % inclination angle error
                obj.N_VAR_GAM_A = varargin{i+1};   
                obj.B_VAR_GAM_A = nan(obj.N_VAR_GAM_A,1);
                % weighted average for gamma_a
                obj.avg_gam_a = 1/obj.N_VAR_GAM_A;     
                % weighted average for variance of gamma_a
                obj.var_gam_a = 1/(obj.N_VAR_GAM_A-1); 
%             elseif strcmp(varargin{i},'N_VAR_GAM_M');
                
                obj.N_VAR_GAM_M = varargin{i+1}; 
                obj.B_VAR_GAM_M = nan(obj.N_VAR_GAM_M,1); % weighted average for gamma_m
                obj.avg_gam_m = 1/obj.N_VAR_GAM_M;        % weighted average for variance of gamma_m
                obj.var_gam_m = 1/(obj.N_VAR_GAM_M-1);                 
            elseif strcmp(varargin{i},'sigma_gyr')
                obj.sigma_gyr = varargin{i+1};
            elseif strcmp(varargin{i},'sigma_acc')
                obj.sigma_acc = varargin{i+1};                
%                 obj.R_gamma_a_lower = ...
%                     varianceInAngleError(obj.sigma_acc,[0 0 obj.gRef]); % determine the lower bound
            else error('Invalid argument');
            end
        end;
        if ~isnan(obj.sigma_gyr),
            % assume noise in the gyro is equal on the x/y axis
            obj.Q = 2*obj.dt^(2)*obj.sigma_gyr^(2); 
            obj.P_pos = obj.Q; % initialise a priori as Q
            obj.P_pri = obj.Q; % initialise a posteriori as Q

            obj.Qmag = 2*obj.dt^(2)*obj.sigma_gyr^(2); 
            obj.Pmag_pos = obj.Qmag; % initialise a priori as Q
            obj.Pmag_pri = obj.Qmag; % initialise a posteriori as Q            
            
            if isnan(obj.N_VAR_GAM_A), error('N_VAR_GAM_A unassigned'); end
%             if isnan(obj.N_VAR_GAM_M), error('N_VAR_GAM_M unassigned'); end
        end
        % --- original assignment of fixed gain
        obj.muAcc  = obj.muAcc*obj.dt; 
        obj.muMag  = obj.muMag*obj.dt;
    end
    %% UpdateNoCallsArcTan
    % 128 multiplications, 43 additions, 43 subtractions, 6 calls to
    % invSqrt (4 multiplications, 2 subtractions)
    % 2 calls to arctan_approx (multiply (x4), addition (x1), subtraction (x2))
    function obj=UpdateNoCallsArcTan(obj,Gyr,Acc,Mag)
        Acc1=Acc(1);    Acc2=Acc(2);    Acc3=Acc(3);
        Gyr1=Gyr(1);    Gyr2=Gyr(2);    Gyr3=Gyr(3);
        Mag1=Mag(1);    Mag2=Mag(2);    Mag3=Mag(3);
        
        muA2 = obj.muAcc/2; 
        muM2 = obj.muMag/2; 
        qG1 = obj.qGlobal(1); 
        qG2 = obj.qGlobal(2); 
        qG3 = obj.qGlobal(3); 
        qG4 = obj.qGlobal(4);            
        %% Correct for gyros - get quaternion from gyro in device frame and rotate
        if ~isnan(Gyr)
            dt2 = 0.5*obj.dt;
            qW1 = 1; qW2 = Gyr1*dt2; qW3 = Gyr2*dt2; qW4 = Gyr3*dt2;
            qWnorm = 1/java.lang.Math.sqrt(qW1*qW1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % can be replaced by the fast inverse square root
            qW1 = qW1*qWnorm; qW2 = qW2*qWnorm; qW3 = qW3*qWnorm; qW4 = qW4*qWnorm;  
            obj.qAngVelDev = [qW1 qW2 qW3 qW4];
            % Convert to back to global frame
            qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4;
            qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3;
            qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2;
            qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; 
        else
            qGi1 = qG1; qGi2 = qG2; qGi3 = qG3; qGi4 = qG4; 
        end
        %% Correct for accelerometer - get next rotation required to align 
        % accelerometer/gravity with 'up'.
        if ~isnan(Acc)
            % intermediate calculation in calculating v' = qvq*
            qa1  = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3;
            qa2  =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2;
            qa3  =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1;
            qa4  =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1;  
            % Convert acceleration to global frame
            gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3;
            gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2;
            gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1;
            obj.gAcc = [gAcc1,gAcc2,gAcc3];
            % Get fraction of rotation from acc in global to up
            if muA2<0,    
                muA2 = 0;    
                warning('Capping mu at 0');
            elseif muA2>0.5,
                muA2 = 0.5;    
                warning('Capping mu at 0.5');
            end
            gAcc2Sq = gAcc2*gAcc2;
            gAcc1Sq = gAcc1*gAcc1;
            gAcc12SumSq = gAcc1Sq+gAcc2Sq;
            axisVecNorm = 1/java.lang.Math.sqrt((gAcc12SumSq)); % can be replaced by the fast inverse square root
            if axisVecNorm == 0
                qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
            else
                gamma_a        = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3);
                obj.gammaA     = gamma_a;
                mu_a_gamma_a   = muA2*gamma_a;
                
                cos_mu_a_gamma_a = 1;            % small angle approximation
                sin_mu_a_gamma_a = mu_a_gamma_a; % small angle approximation
                obj.cos_muA_gammaA = cos_mu_a_gamma_a;
                obj.sin_muA_gammaA = sin_mu_a_gamma_a;
                
                qUp1    = cos_mu_a_gamma_a;
                qUp2    =  gAcc2 * axisVecNorm * sin_mu_a_gamma_a;
                qUp3    = -gAcc1 * axisVecNorm * sin_mu_a_gamma_a; 
                qUpNorm = 1/java.lang.Math.sqrt((qUp1*qUp1 + qUp2*qUp2 + qUp3*qUp3)); % can be replaced by the fast inverse square root
                % normalise to unit quaternions
                qUp1 = qUp1*qUpNorm; qUp2 = qUp2*qUpNorm; qUp3 = qUp3*qUpNorm; qUp4 = 0;
            end
            % Rotate global frame towards 'up'
            qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 - qUp4*qGi4;
            qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 - qUp4*qGi3;
            qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 + qUp4*qGi2;
            qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 + qUp4*qGi1;
        else
            qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
        end
        %% Correct for magnetometer - get next rotation around vertical to align measured global xy 
        % Transform magnetometer reading into global frame
        if ~isnan(Mag)
            qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3;
            qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2;
            qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1;
            qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1;  
            gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3;
            gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2;
            % Get fraction of rotation from mag in global to north (xy components only)
            if muM2<0,    
                muM2 = 0;    warning('Capping mu at 0');
            elseif muM2>0.5,
                muM2 = 0.5;  warning('Capping mu at 0.5');
            end
            qn2 = 0; qn3 = 0;
            if gMag2 == 0 % y-component is 0 therefore pointing north
                qn1 = 1;  qn4 = 0;
            else
                gamma_m      = arctan2_approx(abs(gMag2),gMag1); % absolute value unsure if always correct
                obj.gammaM   = gamma_m;
                mu_m_gamma_m = muM2*gamma_m;

                cos_mu_m_gamma_m = 1;
                sin_mu_m_gamma_m = mu_m_gamma_m;
                
                obj.cos_muM_gammaM = cos_mu_m_gamma_m;
                obj.sin_muM_gammaM = sin_mu_m_gamma_m;
                
                qn1 = cos_mu_m_gamma_m;
                qn4 = java.lang.Math.signum(-gMag2)*sin_mu_m_gamma_m;
                qNorm = 1/java.lang.Math.sqrt(qn1*qn1+qn4*qn4);  % can be replaced by the fast inverse square root
                qn1 = qn1*qNorm;
                qn4 = qn4*qNorm; % normalise to unit quaternions
            end
            % Rotate global frame towards 'north'        
            qGiii1 = qn1*qGii1 - qn2*qGii2 - qn3*qGii3 - qn4*qGii4;
            qGiii2 = qn1*qGii2 + qn2*qGii1 + qn3*qGii4 - qn4*qGii3;
            qGiii3 = qn1*qGii3 - qn2*qGii4 + qn3*qGii1 + qn4*qGii2;
            qGiii4 = qn1*qGii4 + qn2*qGii3 - qn3*qGii2 + qn4*qGii1;
        else
            qGiii1 = qGii1; qGiii2 = qGii2; qGiii3 = qGii3; qGiii4 = qGii4;
        end
        obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];
    end % END UpdateNoCallsArcTan
    
    function obj=UpdateExtCompensation(obj,Gyr,Acc,Mag)
        Acc1=Acc(1);    Acc2=Acc(2);    Acc3=Acc(3);
        Gyr1=Gyr(1);    Gyr2=Gyr(2);    Gyr3=Gyr(3);
        Mag1=Mag(1);    Mag2=Mag(2);    Mag3=Mag(3);
        
        muA2 = obj.muAcc/2; 
        muM2 = obj.muMag/2; 
        qG1 = obj.qGlobal(1); 
        qG2 = obj.qGlobal(2); 
        qG3 = obj.qGlobal(3); 
        qG4 = obj.qGlobal(4);            
        %% Correct for gyros - get quaternion from gyro in device frame and rotate
        if ~isnan(Gyr)
            dt2 = 0.5*obj.dt;
            qW1 = 1; qW2 = Gyr1*dt2; qW3 = Gyr2*dt2; qW4 = Gyr3*dt2;
            qWnorm = 1/java.lang.Math.sqrt(qW1*qW1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % can be replaced by the fast inverse square root
            qW1 = qW1*qWnorm; qW2 = qW2*qWnorm; qW3 = qW3*qWnorm; qW4 = qW4*qWnorm;  
            obj.qAngVelDev = [qW1 qW2 qW3 qW4];
            % Convert to back to global frame
            qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4;
            qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3;
            qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2;
            qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; 
        else
            qGi1 = qG1; qGi2 = qG2; qGi3 = qG3; qGi4 = qG4; 
        end
        %% Correct for accelerometer - get next rotation required to align 
        % accelerometer/gravity with 'up'.
        if ~isnan(Acc)
            % intermediate calculation in calculating v' = qvq*
            qa1  = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3;
            qa2  =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2;
            qa3  =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1;
            qa4  =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1;  
            % Convert acceleration to global frame
            gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3;
            gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2;
            gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1;
            obj.gAcc = [gAcc1,gAcc2,gAcc3];
            % subtract low pass filtered acceleration
            obj.gAcc = [gAcc1,gAcc2,gAcc3] - obj.ext_a.*obj.c_a;
            gAcc1 = obj.gAcc(1);
            gAcc2 = obj.gAcc(2);
            gAcc3 = obj.gAcc(3);
            % Get fraction of rotation from acc in global to up
            if muA2<0,    
                muA2 = 0;    
                warning('Capping mu at 0');
            elseif muA2>0.5,
                muA2 = 0.5;    
                warning('Capping mu at 0.5');
            end
            gAcc2Sq = gAcc2*gAcc2;
            gAcc1Sq = gAcc1*gAcc1;
            gAcc12SumSq = gAcc1Sq+gAcc2Sq;
            axisVecNorm = 1/java.lang.Math.sqrt((gAcc12SumSq)); % can be replaced by the fast inverse square root
            if axisVecNorm == 0
                qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
            else
                gamma_a        = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3);
                obj.gammaA     = gamma_a;
                mu_a_gamma_a   = muA2*gamma_a;
                
                cos_mu_a_gamma_a = 1;            % small angle approximation
                sin_mu_a_gamma_a = mu_a_gamma_a; % small angle approximation
                obj.cos_muA_gammaA = cos_mu_a_gamma_a;
                obj.sin_muA_gammaA = sin_mu_a_gamma_a;
                
                qUp1    = cos_mu_a_gamma_a;
                qUp2    =  gAcc2 * axisVecNorm * sin_mu_a_gamma_a;
                qUp3    = -gAcc1 * axisVecNorm * sin_mu_a_gamma_a; 
                qUpNorm = 1/java.lang.Math.sqrt((qUp1*qUp1 + qUp2*qUp2 + qUp3*qUp3)); % can be replaced by the fast inverse square root
                % normalise to unit quaternions
                qUp1 = qUp1*qUpNorm; qUp2 = qUp2*qUpNorm; qUp3 = qUp3*qUpNorm; qUp4 = 0;
            end
            % Rotate global frame towards 'up'
            qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 - qUp4*qGi4;
            qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 - qUp4*qGi3;
            qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 + qUp4*qGi2;
            qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 + qUp4*qGi1;
        else
            qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
        end
            qGii = [qGii1 qGii2 qGii3 qGii4];
            accGfr = quatmultiply(quatmultiply(qGii, [0 Acc]), quatinv(qGii));
            obj.ext_a = accGfr(2:4) - [0 0 obj.gRef]; % external acceleration in sensor frame
        %% Correct for magnetometer - get next rotation around vertical to align measured global xy 
        % Transform magnetometer reading into global frame
        if ~isnan(Mag)
            qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3;
            qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2;
            qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1;
            qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1;  
            gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3;
            gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2;
            % Get fraction of rotation from mag in global to north (xy components only)
            if muM2<0,    
                muM2 = 0;    warning('Capping mu at 0');
            elseif muM2>0.5,
                muM2 = 0.5;  warning('Capping mu at 0.5');
            end
            qn2 = 0; qn3 = 0;
            if gMag2 == 0 % y-component is 0 therefore pointing north
                qn1 = 1;  qn4 = 0;
            else
                gamma_m      = arctan2_approx(abs(gMag2),gMag1); % absolute value unsure if always correct
                obj.gammaM   = gamma_m;
                mu_m_gamma_m = muM2*gamma_m;

                cos_mu_m_gamma_m = 1;
                sin_mu_m_gamma_m = mu_m_gamma_m;
                
                obj.cos_muM_gammaM = cos_mu_m_gamma_m;
                obj.sin_muM_gammaM = sin_mu_m_gamma_m;
                
                qn1 = cos_mu_m_gamma_m;
                qn4 = java.lang.Math.signum(-gMag2)*sin_mu_m_gamma_m;
                qNorm = 1/java.lang.Math.sqrt(qn1*qn1+qn4*qn4);  % can be replaced by the fast inverse square root
                qn1 = qn1*qNorm;
                qn4 = qn4*qNorm; % normalise to unit quaternions
            end
            % Rotate global frame towards 'north'        
            qGiii1 = qn1*qGii1 - qn2*qGii2 - qn3*qGii3 - qn4*qGii4;
            qGiii2 = qn1*qGii2 + qn2*qGii1 + qn3*qGii4 - qn4*qGii3;
            qGiii3 = qn1*qGii3 - qn2*qGii4 + qn3*qGii1 + qn4*qGii2;
            qGiii4 = qn1*qGii4 + qn2*qGii3 - qn3*qGii2 + qn4*qGii1;
        else
            qGiii1 = qGii1; qGiii2 = qGii2; qGiii3 = qGii3; qGiii4 = qGii4;
        end
        obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];
    end % END UpdateNoCallsArcTan    
    
%% --- CAHRS with indirect Kalman filter for adaptive correction rate
% initial validation conducted by Phils Ngo
%     function obj = Update_IKF(obj,Gyr,Acc,Mag)
%         obj.ctr = obj.ctr + 1;
%         frame = obj.ctr;
%         muM2 = obj.muMag/2;        
%         % FrameNum should be replaced with an internal ctr starting at 1 
%         Acc1=Acc(1);    Acc2=Acc(2);    Acc3=Acc(3);
%         Gyr1=Gyr(1);    Gyr2=Gyr(2);    Gyr3=Gyr(3);
%         Mag1=Mag(1);    Mag2=Mag(2);    Mag3=Mag(3);
% 
%         qG1 = obj.qGlobal(1); 
%         qG2 = obj.qGlobal(2); 
%         qG3 = obj.qGlobal(3); 
%         qG4 = obj.qGlobal(4); 
%         % Correct for gyros - get quaternion from gyro in device frame and rotate
%         if ~isnan(Gyr)
%             dt2 = 0.5*obj.dt;
%             qW1 = 1; qW2 = Gyr1*dt2; qW3 = Gyr2*dt2; qW4 = Gyr3*dt2;
%             qWnorm = 1/java.lang.Math.sqrt(qW1*qW1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % can be replaced by the fast inverse square root
%             qW1 = qW1*qWnorm; qW2 = qW2*qWnorm; qW3 = qW3*qWnorm; qW4 = qW4*qWnorm; 
%             obj.qAngVelDev = [qW1 qW2 qW3 qW4];
%             % Convert to back to global frame
%             qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4;
%             qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3;
%             qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2;
%             qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; 
%         else
%             qGi1 = qG1; qGi2 = qG2; qGi3 = qG3; qGi4 = qG4; 
%         end
%         % Correct for accelerometer - get next rotation required to align 
%         % accelerometer/gravity with 'up'.
%         if ~isnan(Acc)
%             obj.P_pri = obj.P_pos + obj.Q; % 1 addition
%             % intermediate calculation in calculating v' = qvq*
%             qa1  = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3;
%             qa2  =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2;
%             qa3  =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1;
%             qa4  =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1;  
%             % Convert acceleration to global frame
%             gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3;
%             gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2;
%             gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1;
%             obj.gAcc = [gAcc1,gAcc2,gAcc3];
%             % Get fraction of rotation from acc in global to up
%             gAcc2Sq = gAcc2*gAcc2;
%             gAcc1Sq = gAcc1*gAcc1;
%             gAcc12SumSq = gAcc1Sq+gAcc2Sq;
%             axisVecNorm = 1/java.lang.Math.sqrt((gAcc12SumSq)); % can be replaced by the fast inverse square root
%             if axisVecNorm == 0
%                 qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
%             else
%                 gam_a      = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3);
%                 obj.gammaA = gam_a;                
%                 % --- Run time calculation of variance in gamma_a which is
%                 % used to adaptively calculate R_gamma_a
%                 % 1 addition, 1 subtraction
%                 idx      = mod(frame-1, obj.N_VAR_GAM_A)+1; % tells us the current place to input new value and also idx to get rid of
%                 old_gamA = obj.B_VAR_GAM_A(idx); % keep track of oldest value in the B_VAR_GAM_A
%                 if frame > obj.N_VAR_GAM_A
%                     % 3 additions, 4 subtractions, 3 multiplications
%                     obj.new_avg_gam_a = obj.pre_avg_gam_a + obj.avg_gam_a*gam_a - obj.avg_gam_a * old_gamA ;       
%                     obj.var_sum_gam_a = obj.var_sum_gam_a + (gam_a + old_gamA - obj.pre_avg_gam_a - obj.new_avg_gam_a) * (gam_a - old_gamA);
%                 else % the buffer B_VAR_GAM_A is not full yet, use the fixed value
%                     obj.new_avg_gam_a = obj.pre_avg_gam_a * ((frame-1)/frame) + (1/frame)*gam_a;
%                     obj.var_sum_gam_a = obj.var_sum_gam_a + (gam_a - obj.pre_avg_gam_a) * (gam_a - obj.new_avg_gam_a);
%                 end
%                 obj.pre_avg_gam_a   = obj.new_avg_gam_a;
%                 obj.avg_gamma_a     = obj.new_avg_gam_a;
%                 
%                 obj.B_VAR_GAM_A(idx)= gam_a;   % replace oldest idx in B_VAR_GAM_A with newest value
%                 % 1 multiplication
%                 obj.R_gamma_a    = obj.var_sum_gam_a * obj.var_gam_a; % R_gamma_a from Variance of Gamma in Window                  
%                 % 2 addition, 2 multiplication, 1 inverse square root 
% %                 obj.mu_a         = obj.P_pri/(obj.P_pri + obj.R_gamma_a);
%                 S_a_inv          = 1/sqrt( (obj.P_pri + obj.R_gamma_a)*(obj.P_pri + obj.R_gamma_a) );
%                 obj.mu_a         = obj.P_pri*S_a_inv;
%                 % 2 multiplications
%                 mu_a_gamma_a     = (obj.mu_a*gam_a)*0.5; % K*gam_a is the updated error state estimate 
%                 obj.cos_K_gammaA = 1;          % small angle approximation, cos(theta)=1
%                 obj.sin_K_gammaA = mu_a_gamma_a; % small angle approximation, sin(theta)=theta
%                 qUp1             =  obj.cos_K_gammaA;
%                 qUp2             =  gAcc2 * axisVecNorm * obj.sin_K_gammaA;
%                 qUp3             = -gAcc1 * axisVecNorm * obj.sin_K_gammaA; 
%                 qUpNorm          = 1/java.lang.Math.sqrt((qUp1*qUp1 + qUp2*qUp2 + qUp3*qUp3)); % can be replaced by the fast inverse square root
%                 % normalise to unit quaternions
%                 qUp1 = qUp1*qUpNorm; qUp2 = qUp2*qUpNorm; qUp3 = qUp3*qUpNorm; qUp4 = 0;
%                 % update accumulating error
%                 % 1 subtraction, 1 multiplication
%                 obj.P_pos = (1 - obj.mu_a) * obj.P_pri;
%             end
%             % Rotate global frame towards 'up'
%             qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 - qUp4*qGi4;
%             qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 - qUp4*qGi3;
%             qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 + qUp4*qGi2;
%             qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 + qUp4*qGi1;
%         else
%             qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
%         end
%         % Correct for magnetometer - get next rotation around vertical to align measured global xy 
%         % Transform magnetometer reading into global frame
% %         if false%~isnan(Mag)
% %             obj.Pmag_pri = obj.Pmag_pos + obj.Qmag; %% one new plus operation
% %             % intermediate calculation in calculating v' = qvq*
% %             qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3;
% %             qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2;
% %             qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1;
% %             qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1;  
% %             gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3;
% %             gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2;
% % 
% %             qn2 = 0; qn3 = 0;
% %             if gMag2 == 0 % y-component is 0 therefore pointing north
% %                 qn1 = 1;  qn4 = 0;
% %             else
% %                 gamma_m      = arctan2_approx(abs(gMag2),gMag1); 
% %                 obj.gammaM   = gamma_m;
% %                 % --- Run time calculation of variance in gamma_a which is
% %                 % used to adaptively calculate R_gamma_a
% %                 idx      = mod(frame-1, obj.N_VAR_GAM_M)+1; % tells us the current place to input new value and also idx to get rid of
% %                 old_gamM = obj.B_VAR_GAM_M(idx); % keep track of oldest value in the B_VAR_GAM_A
% %                 if frame > obj.N_VAR_GAM_M
% %                     obj.new_avg_gam_m = obj.pre_avg_gam_m + obj.avg_gam_m*gamma_m - obj.avg_gam_m * old_gamM ;       
% %                     obj.var_sum_gam_m = obj.var_sum_gam_m + (gamma_m + old_gamM - obj.pre_avg_gam_m - obj.new_avg_gam_m) * (gamma_m - old_gamM);
% %                 else % the buffer B_VAR_GAM_A is not full yet
% %                     obj.new_avg_gam_m = obj.pre_avg_gam_m * ((frame-1)/frame) + (1/frame)*gamma_m;
% %                     obj.var_sum_gam_m = obj.var_sum_gam_m + (gamma_m - obj.pre_avg_gam_m) * (gamma_m - obj.new_avg_gam_m);
% %                 end
% %                 obj.pre_avg_gam_m   = obj.new_avg_gam_m;
% %                 obj.avg_gamma_m     = obj.new_avg_gam_m;
% %                 
% %                 obj.B_VAR_GAM_M(idx)= gamma_m;   % replace oldest idx in B_VAR_GAM_A with newest value
% %                 obj.R_gamma_m    = obj.var_sum_gam_m * obj.var_gam_m; % R_gamma_a from Variance of Gamma in Window                  
% %                 obj.mu_m         = obj.Pmag_pri/(obj.Pmag_pri + obj.R_gamma_m);                
% %                 mu_m_gamma_m     = (obj.mu_m*gamma_m)/2; % K*gam_a is the updated error state estimate                 
% %                 obj.cos_muM_gammaM = 1;
% %                 obj.sin_muM_gammaM = mu_m_gamma_m;
% %                 qn1 = obj.cos_muM_gammaM;
% %                 qn4 = java.lang.Math.signum(-gMag2)*obj.sin_muM_gammaM;
% %                 qNorm = 1/java.lang.Math.sqrt(qn1*qn1+qn4*qn4);  % can be replaced by the fast inverse square root
% %                 qn1 = qn1*qNorm;
% %                 qn4 = qn4*qNorm; % normalise to unit quaternions
% %             end
% %             % Rotate global frame towards 'north'        
% %             qGiii1 = qn1*qGii1 - qn2*qGii2 - qn3*qGii3 - qn4*qGii4;
% %             qGiii2 = qn1*qGii2 + qn2*qGii1 + qn3*qGii4 - qn4*qGii3;
% %             qGiii3 = qn1*qGii3 - qn2*qGii4 + qn3*qGii1 + qn4*qGii2;
% %             qGiii4 = qn1*qGii4 + qn2*qGii3 - qn3*qGii2 + qn4*qGii1;
% %         else
%             qGiii1 = qGii1; qGiii2 = qGii2; qGiii3 = qGii3; qGiii4 = qGii4;
% %         end
%         obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];        
%     end
    
    function obj = Update_IKF(obj,Gyr,Acc,Mag)
        muM2 = obj.muMag/2;        
        % FrameNum should be replaced with an internal ctr starting at 1 
        Acc1=Acc(1);    Acc2=Acc(2);    Acc3=Acc(3);
        Gyr1=Gyr(1);    Gyr2=Gyr(2);    Gyr3=Gyr(3);
        Mag1=Mag(1);    Mag2=Mag(2);    Mag3=Mag(3);

        qG1 = obj.qGlobal(1); 
        qG2 = obj.qGlobal(2); 
        qG3 = obj.qGlobal(3); 
        qG4 = obj.qGlobal(4); 
        % Correct for gyros - get quaternion from gyro in device frame and rotate
        if ~isnan(Gyr)
            dt2 = 0.5*obj.dt;
            qW1 = 1; qW2 = Gyr1*dt2; qW3 = Gyr2*dt2; qW4 = Gyr3*dt2;
            qWnorm = 1/java.lang.Math.sqrt(qW1*qW1 + qW2*qW2 + qW3*qW3 + qW4*qW4); % can be replaced by the fast inverse square root
            qW1 = qW1*qWnorm; qW2 = qW2*qWnorm; qW3 = qW3*qWnorm; qW4 = qW4*qWnorm; 
            obj.qAngVelDev = [qW1 qW2 qW3 qW4];
            % Convert to back to global frame
            qGi1 = qG1*qW1 - qG2*qW2 - qG3*qW3 - qG4*qW4;
            qGi2 = qG1*qW2 + qG2*qW1 + qG3*qW4 - qG4*qW3;
            qGi3 = qG1*qW3 - qG2*qW4 + qG3*qW1 + qG4*qW2;
            qGi4 = qG1*qW4 + qG2*qW3 - qG3*qW2 + qG4*qW1; 
        else
            qGi1 = qG1; qGi2 = qG2; qGi3 = qG3; qGi4 = qG4; 
        end
        % Correct for accelerometer - get next rotation required to align 
        % accelerometer/gravity with 'up'.
        if ~isnan(Acc)
            % if stationary correct completely,
            
            % if moving adaptive gain with measurement estimation
            % intermediate calculation in calculating v' = qvq*
            qa1  = -qGi2*Acc1 - qGi3*Acc2 - qGi4*Acc3;
            qa2  =  qGi1*Acc1 + qGi3*Acc3 - qGi4*Acc2;
            qa3  =  qGi1*Acc2 - qGi2*Acc3 + qGi4*Acc1;
            qa4  =  qGi1*Acc3 + qGi2*Acc2 - qGi3*Acc1;  
            % Convert acceleration to global frame
            gAcc1 = -qa1*qGi2 + qa2*qGi1 - qa3*qGi4 + qa4*qGi3;
            gAcc2 = -qa1*qGi3 + qa2*qGi4 + qa3*qGi1 - qa4*qGi2;
            gAcc3 = -qa1*qGi4 - qa2*qGi3 + qa3*qGi2 + qa4*qGi1;
            obj.gAcc = [gAcc1,gAcc2,gAcc3];
            % Get fraction of rotation from acc in global to up
            gAcc2Sq = gAcc2*gAcc2;
            gAcc1Sq = gAcc1*gAcc1;
            gAcc12SumSq = gAcc1Sq+gAcc2Sq;
            axisVecNorm = 1/java.lang.Math.sqrt((gAcc12SumSq)); % can be replaced by the fast inverse square root
            if axisVecNorm == 0
                qUp1 = 1; qUp2 = 0; qUp3 = 0; qUp4 = 0;
            else
                gam_a      = arctan2_approx(axisVecNorm*gAcc12SumSq,gAcc3);
                obj.gammaA = gam_a;     
            % Estimate measurement error, gamma_a using a sliding window
% -------------------------------------------------------------------------                
% --- Run time calculation of variance in gamma_a which is
                % used to adaptively calculate R_gamma_a
                % 1 addition, 1 subtraction
                obj.ctr = obj.ctr + 1;
                frame = obj.ctr;
                idx      = mod(frame-1, obj.N_VAR_GAM_A)+1; % tells us the current place to input new value and also idx to get rid of
                old_gamA = obj.B_VAR_GAM_A(idx); % keep track of oldest value in the B_VAR_GAM_A
                if frame > obj.N_VAR_GAM_A
                    % 3 additions, 4 subtractions, 3 multiplications
                    obj.new_avg_gam_a = obj.pre_avg_gam_a + obj.avg_gam_a*gam_a - obj.avg_gam_a * old_gamA ;       
                    obj.var_sum_gam_a = obj.var_sum_gam_a + (gam_a + old_gamA - obj.pre_avg_gam_a - obj.new_avg_gam_a) * (gam_a - old_gamA);
                else % the buffer B_VAR_GAM_A is not full yet, use the fixed value
                    obj.new_avg_gam_a = obj.pre_avg_gam_a * ((frame-1)/frame) + (1/frame)*gam_a;
                    obj.var_sum_gam_a = obj.var_sum_gam_a + (gam_a - obj.pre_avg_gam_a) * (gam_a - obj.new_avg_gam_a);
                end
                obj.pre_avg_gam_a   = obj.new_avg_gam_a;
                obj.avg_gamma_a     = obj.new_avg_gam_a;
                obj.B_VAR_GAM_A(idx)= gam_a;   % replace oldest idx in B_VAR_GAM_A with newest value
                % 1 multiplication
                obj.R_gamma_a    = obj.var_sum_gam_a * obj.var_gam_a; % R_gamma_a from Variance of Gamma in Window                  
                obj.P_pri = obj.P_pos + obj.Q; % 1 addition
                % 2 addition, 2 multiplication, 1 inverse square root 
%                 obj.mu_a         = obj.P_pri/(obj.P_pri + obj.R_gamma_a);
                S_a_inv          = 1/sqrt( (obj.P_pri + obj.R_gamma_a)*(obj.P_pri + obj.R_gamma_a) ); % can be replaced by the fast inverse square root
                obj.mu_a         = obj.P_pri*S_a_inv;
                obj.P_pos        = (1 - obj.mu_a) * obj.P_pri;
                % 2 multiplications
% -------------------------------------------------------------------------                
% or use a fixed gain approach
                
                mu_a_gamma_a     = (obj.mu_a*gam_a)*0.5; % K*gam_a is the updated error state estimate 
                obj.cos_K_gammaA = 1;          % small angle approximation, cos(theta)=1
                obj.sin_K_gammaA = mu_a_gamma_a; % small angle approximation, sin(theta)=theta
                qUp1             =  obj.cos_K_gammaA;
                qUp2             =  gAcc2 * axisVecNorm * obj.sin_K_gammaA;
                qUp3             = -gAcc1 * axisVecNorm * obj.sin_K_gammaA; 
                qUpNorm          = 1/java.lang.Math.sqrt((qUp1*qUp1 + qUp2*qUp2 + qUp3*qUp3)); % can be replaced by the fast inverse square root
                % normalise to unit quaternions
                qUp1 = qUp1*qUpNorm; qUp2 = qUp2*qUpNorm; qUp3 = qUp3*qUpNorm; qUp4 = 0;
                % update accumulating error
                % 1 subtraction, 1 multiplication
                
            end
            % Rotate global frame towards 'up'
            qGii1 = qUp1*qGi1 - qUp2*qGi2 - qUp3*qGi3 - qUp4*qGi4;
            qGii2 = qUp1*qGi2 + qUp2*qGi1 + qUp3*qGi4 - qUp4*qGi3;
            qGii3 = qUp1*qGi3 - qUp2*qGi4 + qUp3*qGi1 + qUp4*qGi2;
            qGii4 = qUp1*qGi4 + qUp2*qGi3 - qUp3*qGi2 + qUp4*qGi1;
        else
            qGii1 = qGi1;  qGii2 = qGi2;  qGii3 = qGi3;  qGii4 = qGi4;
        end
        % Correct for magnetometer - get next rotation around vertical to align measured global xy 
        % Transform magnetometer reading into global frame
        if false%~isnan(Mag)
%             obj.Pmag_pri = obj.Pmag_pos + obj.Qmag; %% one new plus operation
%             % intermediate calculation in calculating v' = qvq*
%             qm1  = -qGii2*Mag1 - qGii3*Mag2 - qGii4*Mag3;
%             qm2  =  qGii1*Mag1 + qGii3*Mag3 - qGii4*Mag2;
%             qm3  =  qGii1*Mag2 - qGii2*Mag3 + qGii4*Mag1;
%             qm4  =  qGii1*Mag3 + qGii2*Mag2 - qGii3*Mag1;  
%             gMag1 = -qm1*qGii2 + qm2*qGii1 - qm3*qGii4 + qm4*qGii3;
%             gMag2 = -qm1*qGii3 + qm2*qGii4 + qm3*qGii1 - qm4*qGii2;
% 
%             qn2 = 0; qn3 = 0;
%             if gMag2 == 0 % y-component is 0 therefore pointing north
%                 qn1 = 1;  qn4 = 0;
%             else
%                 gamma_m      = arctan2_approx(abs(gMag2),gMag1); 
%                 obj.gammaM   = gamma_m;
%                 % --- Run time calculation of variance in gamma_a which is
%                 % used to adaptively calculate R_gamma_a
%                 idx      = mod(frame-1, obj.N_VAR_GAM_M)+1; % tells us the current place to input new value and also idx to get rid of
%                 old_gamM = obj.B_VAR_GAM_M(idx); % keep track of oldest value in the B_VAR_GAM_A
%                 if frame > obj.N_VAR_GAM_M
%                     obj.new_avg_gam_m = obj.pre_avg_gam_m + obj.avg_gam_m*gamma_m - obj.avg_gam_m * old_gamM ;       
%                     obj.var_sum_gam_m = obj.var_sum_gam_m + (gamma_m + old_gamM - obj.pre_avg_gam_m - obj.new_avg_gam_m) * (gamma_m - old_gamM);
%                 else % the buffer B_VAR_GAM_A is not full yet
%                     obj.new_avg_gam_m = obj.pre_avg_gam_m * ((frame-1)/frame) + (1/frame)*gamma_m;
%                     obj.var_sum_gam_m = obj.var_sum_gam_m + (gamma_m - obj.pre_avg_gam_m) * (gamma_m - obj.new_avg_gam_m);
%                 end
%                 obj.pre_avg_gam_m   = obj.new_avg_gam_m;
%                 obj.avg_gamma_m     = obj.new_avg_gam_m;
%                 
%                 obj.B_VAR_GAM_M(idx)= gamma_m;   % replace oldest idx in B_VAR_GAM_A with newest value
%                 obj.R_gamma_m    = obj.var_sum_gam_m * obj.var_gam_m; % R_gamma_a from Variance of Gamma in Window                  
%                 obj.mu_m         = obj.Pmag_pri/(obj.Pmag_pri + obj.R_gamma_m);                
%                 mu_m_gamma_m     = (obj.mu_m*gamma_m)/2; % K*gam_a is the updated error state estimate                 
%                 obj.cos_muM_gammaM = 1;
%                 obj.sin_muM_gammaM = mu_m_gamma_m;
%                 qn1 = obj.cos_muM_gammaM;
%                 qn4 = java.lang.Math.signum(-gMag2)*obj.sin_muM_gammaM;
%                 qNorm = 1/java.lang.Math.sqrt(qn1*qn1+qn4*qn4);  % can be replaced by the fast inverse square root
%                 qn1 = qn1*qNorm;
%                 qn4 = qn4*qNorm; % normalise to unit quaternions
%             end
%             % Rotate global frame towards 'north'        
%             qGiii1 = qn1*qGii1 - qn2*qGii2 - qn3*qGii3 - qn4*qGii4;
%             qGiii2 = qn1*qGii2 + qn2*qGii1 + qn3*qGii4 - qn4*qGii3;
%             qGiii3 = qn1*qGii3 - qn2*qGii4 + qn3*qGii1 + qn4*qGii2;
%             qGiii4 = qn1*qGii4 + qn2*qGii3 - qn3*qGii2 + qn4*qGii1;
        else
            qGiii1 = qGii1; qGiii2 = qGii2; qGiii3 = qGii3; qGiii4 = qGii4;
        end
        obj.qGlobal = [qGiii1 qGiii2 qGiii3 qGiii4];        
    end    
    
end 
end

function [atan2x]=arctan2_approx(y,x)
    % only works when y is positive i.e. the numerator in atan2(y,x)
    if x >= 0
        if x>=y
            invsqrtxsq = java.lang.Math.signum(x)*(1/java.lang.Math.sqrt(x*x)); % can be replaced by the fast inverse square root
            %atan2x = arctan_approx(y/x);
            atan2x = arctan_approx(y*invsqrtxsq);
        else % x < y
            invsqrtysq = java.lang.Math.signum(y)*(1/java.lang.Math.sqrt(y*y)); % can be replaced by the fast inverse square root
            %atan2x = pi/2 - arctan_approx(x/y);
            atan2x = pi/2 - arctan_approx(x*invsqrtysq);
        end
    else % x<=0
        if y>abs(x)
            %atan2x = pi/2 + arctan_approx(abs(x)/y);
            invsqrtysq = java.lang.Math.signum(y)*(1/java.lang.Math.sqrt(y*y)); % can be replaced by the fast inverse square root
            atan2x = pi/2 + arctan_approx(abs(x)*invsqrtysq);
        else 
            %atan2x = pi - arctan_approx(y/abs(x));                        
            invsqrtxsq = 1/java.lang.Math.sqrt(x*x); % can be replaced by the fast inverse square root
            atan2x = pi - arctan_approx(y*invsqrtxsq);
        end
    end
end

function atanx = arctan_approx(x)
qtr_pi = pi/4;
% multiply (x4)
% addition (x1)
% subtraction (x2)
atanx = (qtr_pi.*x) - x.*(abs(x)-1) .* (0.2447 + 0.0663 .* abs(x));
end

%
function [varErr] = varianceInAngleError(sigma_acc,refVec)
% sigma_acc - std deviation in the measurements from the accelerometer
% refVec - the reference vector [0 0 1] if accelerometer
rng(1);
N_SIMS = 10000;	
upVec = bsxfun(@plus,refVec,normrnd(0, sigma_acc,N_SIMS,3));
upRef = repmat(refVec, N_SIMS,1);
cos_theta = NaN(N_SIMS,1);
angle_rad = NaN(N_SIMS,1);

for a = 1:N_SIMS
    cos_theta(a) = dot(upVec(a,:),upRef(a,:))...
                        /(norm(upVec(a,:))*norm(upRef(a,:)));
    angle_rad(a) = acos(cos_theta(a));
end
varErr = var(angle_rad);
%%
% updateFigureContents('Histograms');
% subplot(1,4,1)
% [nx,cx]=
% hist(upVec(:,1));
% subplot(1,4,2)
% [ny,cy]=
% hist(upVec(:,2));
% subplot(1,4,3)
% [nz,cz]=
% hist(upVec(:,3));
% subplot(1,4,4)
% [na,ca]=
% hist(angle_rad);
% bar(na,ca,1)
%%
% updateFigureContents('upAngleVariance');hold on;axis('square');grid on;
% surbplot(2,1,1)
% plot3(upVec(:,1),upVec(:,2),upVec(:,3),'x');
% surf(upVec(:,1),upVec(:,2),angle_rad);
% 
% line([0 refVec(1)],[0 refVec(2)],[0 refVec(3)],'Color','k');
% xlabel('x');ylabel('y');zlabel('z');
% plot3(upVec(:,1),upVec(:,2),upVec(:,3),'x');
end
