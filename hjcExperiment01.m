%% Hinge Joint Constraint Experiment 01
% Author: Luke Sy
%
% play with the hinge joint constraint and understand its behavior
% 
% Findings:
% 1. The update ('d') is only based on position error. If P=I, only the
% position states will be updated depending on P (higher value means less 
% trust worthy. The velocity states will only be updated due to covariance 
% between position and velocity.
% 2. The update does not change the knee angle (alpha). Hence, if position
% error causes knee angle error, the error is carried over and not fixed by
% this constraint update.

d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
D = [-eye(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
     -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3)];
 
PELV_CS = eye(3,3);
LTIB_CS = eye(3,3);
RTIB_CS = eye(3,3);
P = 0.1*eye(18,18);
P(1:3,1:3) = 100*eye(3,3);
% P(1:3,4:6) = eye(3,3);
% P(4:6,1:3) = eye(3,3);

% Ideal case straight
% xhat = [0 0 0 0.1 0.1 0 ...
%         0 0.25 -1.1 0.1 -0.1 0 ...
%         0 -0.25 -1.1 0.1 -0.1 0]';
% Non-ideal case straight, longer/shorter length
xhat = [0 0 0 0.1 0.1 0 ...
        0 0.3 -1.2 0.1 -0.1 0 ...
        0 -0.3 -1.0 0.1 -0.1 0]';

% Ideal case bent
% xhat = [0 0 0 0.1 0.1 0 ...
%         0.6 0.25 -0.5 0.1 -0.1 0 ...
%         0.6 -0.25 -0.5 0.1 -0.1 0]';
% Non-ideal case bent, longer/shorter length
% xhat = [0 0 0 0.1 0.1 0 ...
%         0.7 0.25 -0.7 0.1 -0.1 0 ...
%         0.7 -0.25 -0.4 0.1 -0.1 0]';

% Ideal case slightly bent
% xhat = [0 0 0 0.1 0.1 0 ...
%         0.6*sin(pi/4) 0.25 0.6*cos(pi/4)-0.5 0.1 -0.1 0 ...
%         0.6*sin(pi/4) -0.25 0.6*cos(pi/4)-0.5 0.1 -0.1 0]';
% Ideal case slightly bent, knee not hinge joint
% xhat = [0 0 0 0.1 0.1 0 ...
%         0.6*sin(pi/4) 0.5 0.6*cos(pi/4)-0.5 0.1 -0.1 0 ...
%         0.6*sin(pi/4) -0.5 0.6*cos(pi/4)-0.5 0.1 -0.1 0]';

LKNE = xhat(7:9,1) + d_ltibia*LTIB_CS(:,3);
RKNE = xhat(13:15,1) + d_rtibia*RTIB_CS(:,3);
LFEM_z = xhat(1:3,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
RFEM_z = xhat(1:3,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;
alpha_lk = acos(dot(LFEM_z, LTIB_CS(:,3))/(norm(LFEM_z)*norm(LTIB_CS(:,3))));
alpha_rk = acos(dot(RFEM_z, RTIB_CS(:,3))/(norm(RFEM_z)*norm(RTIB_CS(:,3))));

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

LKNE = xhat(7:9,1) + d_ltibia*LTIB_CS(:,3);
RKNE = xhat(13:15,1) + d_rtibia*RTIB_CS(:,3);
LHIP = xhat(1:3,1)+d_pelvis/2*PELV_CS(:,2);
RHIP = xhat(1:3,1)-d_pelvis/2*PELV_CS(:,2);
LFEM_z = LHIP-LKNE;
RFEM_z = RHIP-RKNE;
alpha_lk = acos(dot(LFEM_z, LTIB_CS(:,3))/(norm(LFEM_z)*norm(LTIB_CS(:,3))));
alpha_rk = acos(dot(RFEM_z, RTIB_CS(:,3))/(norm(RFEM_z)*norm(RTIB_CS(:,3))));