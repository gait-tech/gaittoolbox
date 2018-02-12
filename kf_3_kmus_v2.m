function [ varargout ] = kf_3_kmus_v2(fs, ...
    sigma_acc_MP, sigma_acc_LA, sigma_acc_RA, P0, ...
    x0_pos_MP, x0_vel_MP, gfr_acc_MP, bIsStatMP, q_MP, ...
    x0_pos_LA, x0_vel_LA, gfr_acc_LA, bIsStatLA, q_LA, ...
    x0_pos_RA, x0_vel_RA, gfr_acc_RA, bIsStatRA, q_RA, ...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
    zerovel_update, uwb_update, kneecoplanar_constraint, ...
    femurdist_constraint, kneeangle_constraint)
% KF_3_KMUS Kalman Filter for performing sensor fusion on the trajectory of
% three KMUs presumably worn on the body in the following configuration: mid
% pelvis, left ankle, right ankle
% In this state space model, the position and velocity of each kinematic
% measurement unit (KMU) is estimated in 3D space by combining the
% information from each KMU in a kalman filter. NOTE: pay special attention 
% to units:
% position (meters)
% velocity (m/s)
% acceleration (m/2^2)
% uwb_mea (meters)
%
% Author: Michael Del Rosario, Luke Wicent Sy
%
% Inputs::
%   fs - sampling frequency of the magnetic and inertial measurement units
%   sigma_acc - user specified process noise, i.e., the standard deviation
%               in the accelerometer measurements when subjected to a known
%               acceleration
%   x0_pos_MP  - the initial estimate of the mid-pelvis in the GFR
%
%   x0_vel_MP  - the initial velocity of the mid-pelvis in the GFR
%
%   gfr_acc_MP - the acceleration of the mid-pelvis in the GFR
%
%   x0_pos_LA  - the initial estimate of the left ankle in the GFR
%
%   x0_vel_LA  - the initial velocity of the left ankle in the GFR
%
%   gfr_acc_LA - the acceleration of the left ankle in the GFR
%
%   x0_pos_RA  - the initial estimate of the right ankle in the GFR
%
%   x0_vel_RA  - the initial velocity of the right ankle in the GFR
%
%   gfr_acc_RA - the acceleration of the right ankle in the GFR
%
%   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%                
%   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%
%   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   q_MP       - mid  pelvis orientation in the GFR (quaternion)
%   q_LA       - left  ankle orientation in the GFR (quaternion)
%   q_RA       - right ankle orientation in the GFR (quaternion)
%   d_pelvis   - pelvis width
%   d_rfemur   - right femur length
%   d_lfemur   - left femur length
%   d_rtibia   - right tibia length
%   d_ltibia   - left tibia length
%   uwb_mea    - a structure containing the range measurements (m) between 

idx_pos_MP = 1:3; % column idx corresponding to the mid-pelvis position
idx_vel_MP = 4:6; % column idx corresponding to the mid-pelvis velocity
idx_pos_LA = 7:9; % column idx corresponding to the left ankle position
idx_vel_LA = 10:12; % column idx corresponding to the left ankle velocity
idx_pos_RA = 13:15; % column idx corresponding to the right ankle position
idx_vel_RA = 16:18; % column idx corresponding to the right ankle velocity

% pseudo UWB measurements corresponding to the euclidean distance between
% pairs of KMUs.
% NOTE: the column order of these measurements is important. Assume the
% following order unless stated otherwise:
%  uwb_MP_LA_RA =  ['mid-pelvis to left ankle',
%                   'mid-pelvis to right ankle',
%                   'left ankle to right ankle'];

uwb_MP_LA_RA = [uwb_mea.left_tibia_mid_pelvis,...
                uwb_mea.mid_pelvis_right_tibia,...
                uwb_mea.left_tibia_right_tibia];

% initialise state vector (must be column)
x0 = [x0_pos_MP, x0_vel_MP, x0_pos_LA, x0_vel_LA, x0_pos_RA, x0_vel_RA]';
xhat = x0;

% specify the measurement noise in the UWB measurements, these may be
% different for each KMU sensor pair. It is likely that the range between
% feet/ankles will be the most accurate due to less "no line of sight"
% periods. Note: units on sigma_uwb = meters

sigma_uwb_mp_la = 0.2;
sigma_uwb_mp_ra = 0.2;
sigma_uwb_la_ra = 0.1;

R_uwb = zeros(3,3);
R_uwb(1,1) = sigma_uwb_mp_la*sigma_uwb_mp_la;
R_uwb(2,2) = sigma_uwb_mp_ra*sigma_uwb_mp_ra;
R_uwb(3,3) = sigma_uwb_la_ra*sigma_uwb_la_ra;
% if all measurements have the same measurement model, replace with
% R_uwb = eye(3).* sigma_uwb^2;

% specify the noise level in the psuedo zero-velocity measurements, this
% should be equal for both ankles
sigma_vel = 0.5; % (units m/s)

dt = 1/fs;        % assume constant sampling interval
dt2 = 0.5*dt*dt;  % local variable for readability
N_STATES = 18;
I_18 = eye(N_STATES);
% state transition matrix encodes the relationship between previous state
% estimate and current state estimate
A = eye(N_STATES,N_STATES);
% velocity from previous time-step affects position at next time-step, i.e.,
% x = x(t-1) + v(t-1)*dt
A(1:3,4:6)     = dt.*eye(3); % mid pelvis
A(7:9,10:12)   = dt.*eye(3); % left ankle
A(13:15,16:18) = dt.*eye(3); % right ankle
% calculate the control input model, B, i.e., the matrix applied to the
% control/input vector, in this case the accelerometer measurements in 3D
% x = x(t-1) + v(t-1)*dt + 0.5*a(t)*dt^2
%    [  A * xhat       ] + [ B * acc]
B            = zeros(18,9);
B(1:3,1:3)   = dt2.*eye(3);
B(4:6,1:3)   = dt .*eye(3);
B(7:9,4:6)   = dt2.*eye(3);
B(10:12,4:6) = dt .*eye(3);
B(13:15,7:9) = dt2.*eye(3);
B(16:18,7:9) = dt .*eye(3);

% Constraint
D = [-eye(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
     -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3)];

% Initialise process noise covariance
Q = zeros(9,9);
Q(1:3,1:3) = eye(3,3)*sigma_acc_MP^2;
Q(4:6,4:6) = eye(3,3)*sigma_acc_LA^2;
Q(7:9,7:9) = eye(3,3)*sigma_acc_RA^2;
Q = B * Q * B';

% initialise covariance in the state estimate
if islogical(P0) && ~P0
    P = Q;
else
    P = P0;
end

% check that all accelerometer measurements are equal dimensions
[N_MP,~] = size(gfr_acc_MP);
[N_LA,~] = size(gfr_acc_LA);
[N_RA,~] = size(gfr_acc_RA);

if (N_MP ~= N_LA) || (N_MP ~= N_RA) || (N_LA ~= N_RA), 
    error('Accelerometer measurements from all sensors unequal'); 
end

% local variable assignment for readability
acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
% allocate memory to store apriori and aposteriori state estimates, xhat,
% and error covariances in the state estimate, P_pri, P_pos
xhat_pri = nan(N_MP,     N_STATES);
P_pri    = nan(N_STATES, N_STATES, N_MP);

xhat_pos = nan(N_MP,     N_STATES);
P_pos    = nan(N_STATES, N_STATES, N_MP);

tmp_dat = struct;
tmp_dat.LFEO = nan(N_MP, 3); tmp_dat.RFEO = nan(N_MP, 3);
tmp_dat.LFEP = nan(N_MP, 3); tmp_dat.RFEP = nan(N_MP, 3);
tmp_dat.qLFemur = nan(N_MP, 4); tmp_dat.qRFemur = nan(N_MP, 4);
tmp_dat.predState = nan(N_MP, N_STATES);
% if zerovel_update
    tmp_dat.zuptStateL = nan(N_MP, N_STATES);
    tmp_dat.zuptStateR = nan(N_MP, N_STATES);
% end
% if uwb_update
    tmp_dat.uwbuptState = nan(N_MP, N_STATES);
% end
% if kneecoplanar_constraint
    tmp_dat.cpkneeState = nan(N_MP, N_STATES);
    tmp_dat.cpkneeStateKk = nan(N_STATES, 6, N_MP);
    tmp_dat.cpkneeStateRes = nan(N_MP, 6);
% end
if femurdist_constraint
    tmp_dat.femdistState = nan(N_MP, N_STATES);
end

for n = 1:N_MP
%% -----------------------------------------------------------------------
%  ---- Prediction Step using accelerometer measurements ----    
    xhat = A * xhat + B * acc(n,:)' ;
    P_min= A * P * A' + Q;
    xhat_pri(n,:) = xhat;
    P_pri(:,:,n)  = P_min;
    
    tmp_dat.predState(n,:) = xhat;
    
%% ------------------------------------------------------------------------
%  ---- Implement the zero velocity update step
%  matrices beginnning with 'H_' are the 'observation matrices' that map
%  the variables in the state estimate vector, xhat, to the measurement
%  domain. In this case we are using
    if zerovel_update
        ctrZUPT = 0;
        H_MP=[];
        if bIsStatMP(n), 
            ctrZUPT = ctrZUPT+1; 
            H_MP = zeros(3,N_STATES);
            H_MP(:,idx_vel_MP) = eye(3);
        end
        H_LA=[];
        if bIsStatLA(n), 
            ctrZUPT = ctrZUPT+1; 
            H_LA = zeros(3,N_STATES);
            H_LA(:,idx_vel_LA) = eye(3);
        end
        H_RA=[];
        if bIsStatRA(n), 
            ctrZUPT = ctrZUPT+1; 
            H_RA = zeros(3,N_STATES);
            H_RA(:,idx_vel_RA) = eye(3);
        end

        N_PSEUDO_MEA = 3*ctrZUPT; % the number of pseudo measurements
        if ctrZUPT > 0
            H_zup = [H_MP;H_LA;H_RA];
            z_zup = zeros(N_PSEUDO_MEA,1);
            y_zup = z_zup - H_zup*xhat;

            R_zup = eye(N_PSEUDO_MEA).*sigma_vel^2;

            S_zup = H_zup*P_min*H_zup'+ R_zup;
            K_zup = P_min*H_zup' /(S_zup);%P_min*H_zup' * inv(S_zup);

            xhat = xhat + K_zup*y_zup;
            P_min = (I_18 - K_zup*H_zup)*P_min;
            
        end
        
        if bIsStatLA(n)
            tmp_dat.zuptStateL(n,:) = xhat;
        end
        if bIsStatRA(n)
            tmp_dat.zuptStateR(n,:) = xhat;
        end
    end
%% ---- Kalman Filter Update Step using UWB measurements ---- 
% this correction step should be done last
    if uwb_update
        diff_MP_LA = xhat(idx_pos_MP)'-xhat(idx_pos_LA)';
        diff_MP_RA = xhat(idx_pos_MP)'-xhat(idx_pos_RA)';
        diff_LA_RA = xhat(idx_pos_LA)'-xhat(idx_pos_RA)';
        % the observation model
        h_uwb_est = [vecnormalize( diff_MP_LA );
                     vecnormalize( diff_MP_RA );
                     vecnormalize( diff_LA_RA );];
        % calculate the measurement residual, i.e., the difference between the
        % predicted measurements from the observation model, h(xhat), and the
        % measurements obtained directly from the UWB sensors.
        y_innovation = uwb_MP_LA_RA(n,:)' - h_uwb_est;

        H = zeros(3, N_STATES);
        % Jacobian of observation model with respect to elements of the state
        % estimate xhat, the order of these matrix elements MUST be preserved
        H(1,1) = diff_MP_LA(1)/h_uwb_est(1);
        H(1,2) = diff_MP_LA(2)/h_uwb_est(1);
        H(1,3) = diff_MP_LA(3)/h_uwb_est(1);
        H(1,7) = -diff_MP_LA(1)/h_uwb_est(1);
        H(1,8) = -diff_MP_LA(2)/h_uwb_est(1);
        H(1,9) = -diff_MP_LA(3)/h_uwb_est(1);    

        H(2,1)  =  diff_MP_RA(1)/h_uwb_est(2);
        H(2,2)  =  diff_MP_RA(2)/h_uwb_est(2);
        H(2,3)  =  diff_MP_RA(3)/h_uwb_est(2);
        H(2,13) = -diff_MP_RA(1)/h_uwb_est(2);
        H(2,14) = -diff_MP_RA(2)/h_uwb_est(2);
        H(2,15) = -diff_MP_RA(3)/h_uwb_est(2);

        H(3,7)  =  diff_LA_RA(1)/h_uwb_est(3);
        H(3,8)  =  diff_LA_RA(2)/h_uwb_est(3);
        H(3,9)  =  diff_LA_RA(3)/h_uwb_est(3);
        H(3,13) = -diff_LA_RA(1)/h_uwb_est(3);
        H(3,14) = -diff_LA_RA(2)/h_uwb_est(3);
        H(3,15) = -diff_LA_RA(3)/h_uwb_est(3);    
        % Calculate the covariance in the measurement residual
        S    = ((H * P_min) * H') + R_uwb; % scalar
        % Calculate the kalman gain
        K    = P_min * H' * S^(-1);
        % Update state estimate
        xhat = xhat + K * y_innovation;
        % Update Covariance in the state estimate
        P    = (I_18 - K * H) * P_min;
        
        tmp_dat.uwbuptState(n,:) = xhat;
    else
        P = P_min;
    end
    
%% -----------------------------------------------------------------------
%  ---- Constraint update step using anthropometric measurement ----  
    % CS initializations
    LTIB_CS = quat2rotm(q_LA(n,:));
    RTIB_CS = quat2rotm(q_RA(n,:));
    PELV_CS = quat2rotm(q_MP(n,:));
    
    if kneecoplanar_constraint
        % calculate the location of the knee
        LKNE = xhat(7:9,1) + d_ltibia*LTIB_CS(:,3);
        RKNE = xhat(13:15,1) + d_rtibia*RTIB_CS(:,3);

        % calculate the z axis of the femur
        LFEM_z = xhat(1:3,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
        RFEM_z = xhat(1:3,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;

        % calculate the z axis of the tibia
        LTIB_z = LTIB_CS(:,3);
        RTIB_z = RTIB_CS(:,3);

        % calculate alpha_lk and alpha_rk
        alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
        alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));
        
        % setup the constraint equations
        d_k = [ (d_pelvis/2*PELV_CS(:,2) ...
                 -d_lfemur*cos(alpha_lk)*LTIB_CS(:,3) ...
                 +d_lfemur*sin(alpha_lk)*LTIB_CS(:,1) ...
                 -d_ltibia*LTIB_CS(:,3)) ; ...
                (-d_pelvis/2*PELV_CS(:,2)+ ...
                 -d_rfemur*cos(alpha_rk)*RTIB_CS(:,3) ...
                 +d_rfemur*sin(alpha_rk)*RTIB_CS(:,1) ...
                 -d_rtibia*RTIB_CS(:,3)) ];

        res = d_k-D*xhat;
        Kk = P*D'*(D*P*D')^(-1);
        dx = Kk*(res);
        xhat = xhat + dx;
        
        tmp_dat.cpkneeState(n,:) = xhat;
        tmp_dat.cpkneeStateRes(n,:) = res;
        tmp_dat.cpkneeStateKk(:,:,n) = Kk;
    end
    
    if femurdist_constraint
        diff_LH_LK = xhat(idx_pos_MP)+d_pelvis/2*PELV_CS(:,2)...
            -xhat(idx_pos_LA)-d_ltibia*LTIB_CS(:,3);
        diff_RH_RK = xhat(idx_pos_MP)-d_pelvis/2*PELV_CS(:,2)...
            -xhat(idx_pos_RA)-d_rtibia*RTIB_CS(:,3);
        
        % the observation model
        dist_est = [norm(diff_LH_LK);
                    norm(diff_RH_RK)];

        D = zeros(2, N_STATES);
        % Jacobian of observation model with respect to elements of the state
        % estimate xhat, the order of these matrix elements MUST be preserved
        D(1,1) = diff_LH_LK(1)/dist_est(1);
        D(1,2) = diff_LH_LK(2)/dist_est(1);
        D(1,3) = diff_LH_LK(3)/dist_est(1);
        D(1,7) = -diff_LH_LK(1)/dist_est(1);
        D(1,8) = -diff_LH_LK(2)/dist_est(1);
        D(1,9) = -diff_LH_LK(3)/dist_est(1);    

        D(2,1)  =  diff_RH_RK(1)/dist_est(2);
        D(2,2)  =  diff_RH_RK(2)/dist_est(2);
        D(2,3)  =  diff_RH_RK(3)/dist_est(2);
        D(2,13) = -diff_RH_RK(1)/dist_est(2);
        D(2,14) = -diff_RH_RK(2)/dist_est(2);
        D(2,15) = -diff_RH_RK(3)/dist_est(2);
        
        d_k = [d_lfemur; d_rfemur] - dist_est + D*xhat;

        dx = P*D'*(D*P*D')^(-1)*(d_k-D*xhat);
        xhat = xhat + dx;
        
        tmp_dat.femdistState(n,:) = xhat;
    end
    
    if kneeangle_constraint
        
    end
    
    tmp_dat.LFEO(n,:) = xhat(idx_pos_LA)+d_ltibia*LTIB_CS(:,3);
    tmp_dat.RFEO(n,:) = xhat(idx_pos_RA)+d_rtibia*RTIB_CS(:,3);
    tmp_dat.LFEP(n,:) = xhat(idx_pos_MP)+d_pelvis/2*PELV_CS(:,2);
    tmp_dat.RFEP(n,:) = xhat(idx_pos_MP)-d_pelvis/2*PELV_CS(:,2);
    LFEM_z = (tmp_dat.LFEP(n,:)-tmp_dat.LFEO(n,:))';
    LFEM_y = LTIB_CS(:,2);
    LFEM_x = cross(LFEM_y, LFEM_z);
    RFEM_z = (tmp_dat.RFEP(n,:)-tmp_dat.RFEO(n,:))';
    RFEM_y = RTIB_CS(:,2);
    RFEM_x = cross(RFEM_y, RFEM_z);
    tmp_dat.qLTH(n,:) = rotm2quat([LFEM_x LFEM_y LFEM_z]);
    tmp_dat.qRTH(n,:) = rotm2quat([RFEM_x RFEM_y RFEM_z]);
    
    xhat_pos(n,:) = xhat;
    P_pos(:,:,n)  = P;
end


allIdx = 1:N_MP; 
arr_markers = xhat_pri;
vicX = reshape(arr_markers(allIdx,1:3:end),[],1); XLIM = [min(vicX),max(vicX)];
vicY = reshape(arr_markers(allIdx,2:3:end),[],1); YLIM = [min(vicY),max(vicY)];
vicZ = reshape(arr_markers(allIdx,3:3:end),[],1); ZLIM = [min(vicZ),max(vicZ)];

% figure;grid on;
% xlabel('x - Forward (m)');ylabel('y - East (m)');zlabel('z - Vertical (m)');
% hold on;grid on;axis('equal');view([-51 12]);
% plot3(xhat_pos(:,1),xhat_pos(:,2),xhat_pos(:,3),'k');
% plot3(xhat_pos(:,7),xhat_pos(:,8),xhat_pos(:,9),'b');
% plot3(xhat_pos(:,13),xhat_pos(:,14),xhat_pos(:,15),'r');
% plot3(tmp_dat(:,1),tmp_dat(:,2),tmp_dat(:,3),'g');
% plot3(tmp_dat(:,4),tmp_dat(:,5),tmp_dat(:,6),'v');
% xlim(XLIM);ylim(YLIM);zlim(ZLIM);

if nargout >=1
    varargout{1} = xhat_pri;    
end

if nargout >=2
    varargout{2} = xhat_pos;
end

if nargout >=3
    varargout{3} = tmp_dat;
end


end

