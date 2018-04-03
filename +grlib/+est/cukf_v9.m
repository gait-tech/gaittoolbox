%% cukf_v9 (added velocity-level knee constraint

%v10 will include accel bias in state/zero-velocity-update

% ***SPECIAL NOTE FOR v9: x containts [pos, vel, accel, quat, ang. vel.]' of each sensor
%resulting in a state vec of 48 elements (3 sensors, 16 states per sensor)
%order of sensors in state is MP, LA, RA
% z now containts [accel, quat, ang. vel]' x3 in order MP, LA, RA
% therefore, z has length = 3x10 = 30

function [ x_rec, xa_rec, qFEM ] = cukf_v9(x,P,Q,R,N_MP,nMeas,acc,fs,...
    q_MP, q_LA, q_RA, w_MP_gfr__s, w_LA_gfr__s, w_RA_gfr__s, d_pelvis, d_lfemur, d_rfemur,...
    d_ltibia, d_rtibia,isConstr)
%% Unscented Kalman filter for state estimation of human lower limbs (a
% nonlinear dynamical system)
% [x_rec, xa_rec, qFEM] = ukf(x,P,Q,R,N_MP,nMeas,acc,fs,...
%    q_MP, q_LA, q_RA, d_pelvis, d_lfemur, d_rfemur,...
%   d_ltibia, d_rtibia,isConstr) assumes additive noise
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   x: state at time K
%           P: "a priori" estimated state covariance
%           Q: initial estimated process noise covariance
%           R: initial estimated measurement noise covariance           z: current measurement
%           N_MP: number of timesteps to step forward using filter
%           nMeas: number of dimensions in measurment vec (mes. space)
%           acc: acceleration measuremnt matrix
%           fs: sampleing rate from sensors (1/fs = tiemstep length)
%           q__:q_MP,q_LA,q_RA quaternions of midpel, left ankle, right
%               ankle in gfr
%           w_: w_segment_rel2frame__describedInFrame: angular velocity of each sensor relative to gfr, described in sensor frame
%           d__:d_rtib, d_lfib, d_Lfem, d_rfem, d_pel is length of interest
%           (principle axis) in each segment
% Output:   x_rec: "a posteriori" state estimate
%           xa_rec: "a posteriori" augmented state components
%           qFEM: quaternion representing orientation of femur in gfr

L = length(x);                                 %numer of states
m = nMeas;
alpha = 1e-4;     %should be small (0<alpha<1)                            %default, tunable
ki = 0;            %start with 0 change if needed                           %default, tunable
beta = 2;     % assuming gaussian then use 2                                %default, tunable
%hjc = 3;%3 or 4 - see hingeJoint_cosntr func. for details
lambda = alpha^2*(L+ki)-L;                    %constant
c = L+lambda;                                 %constant
Wm = [lambda/c 0.5/c*ones(1,2*L)];           %weight means
Wc = Wm;
Wc(1) = Wc(1)+(1-alpha^2+beta);               %weight covariance
mu = sqrt(c);
S = chol(P)';
%S = eye(L);
%disp(size(S))
x_rec = nan(L,N_MP);
xa_rec = nan(12,N_MP); % 12 = femur proximal and distal points, both legs(2x2x3)
qFEM = nan(8,N_MP); %8 = 4 per quaterinoin and 2 legs
%vanderwerve 2001 paper wm, wc, and constant definitions

for k=1:N_MP
        disp('k')
        disp(k)
    z = [acc(k,1:3)'; q_MP(k,:)'; w_MP_gfr__s(k,:)'; acc(k,4:6)'; q_LA(k,:)'; w_LA_gfr__s(k,:)'; acc(k,7:9)'; q_RA(k,:)'; w_RA_gfr__s(k,:)'];
    %could constrain raw measurements here
    
    X = genSig(x,S,mu);                            %sigma points around x
    
    [x1,X_post,Px,X_dev,S] = ut(@f_pvaqw,X,Wm,Wc,L,Q,fs);          %unscented transformation of process
    
    [z1,Z_post,Pz,Z_dev,Sy] = ut(@h_pvaqw,X_post,Wm,Wc,m,R,fs);       %unscented transformation of measurments
    %constrain measuremnts
    Pxz=X_dev*diag(Wc)*Z_dev';                        %transformed cross-covariance
    
    K=(Pxz/Sy)/Sy';
    %K=Pxz/(Pz);
    %     disp('dim K')
    %     disp(size(K));
    P = Px-K*Pxz';
    % if isConstr
    %         xhat = x1;
    %         options = optimoptions('fmincon','Algorithm','sqp','Display','off',...
    %             'OptimalityTolerance', 1e-4, 'ConstraintTolerance', 1e-4,...
    %             'MaxFunctionEvaluations',15000); % run interior-point algorithm
    %     x1 = fmincon(@(x1) L2Dist(x1,xhat,S),xhat,[],[],[],[],[],[],@(x1) hingeJoint_constrNL_q(x1,...
    %     d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia),options);
    % end
    x=x1+K*(z-z1);                              %measurement update
    
    if isConstr
        xhat = x;
        options = optimoptions('fmincon','Algorithm','sqp','Display','off',...
            'OptimalityTolerance', 1e-3, 'ConstraintTolerance', 1e-3,...
            'MaxFunctionEvaluations',15000); % run interior-point algorithm
        x = fmincon(@(x) L2Dist(x,xhat,S),xhat,[],[],[],[],[],[],@(x) hingeJoint_constrNL_q(x,...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia),options);
        
        %         [x] = hingeJoint_constrNL(x,hjc,P,k,q_MP, q_LA, q_RA,...
        %             d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia);
    end
    U = K*Sy';
    for i = 1:m
        S = cholupdate(S, U(:,i), '-');
    end
    %% gen augmented x return values
    idx_pos_MP = [1:3]';
    idx_pos_LA = [17:19]';
    idx_pos_RA = [33:35]';
    idx_q_MP = [10:13]';
    idx_q_LA = [26:29]';
    idx_q_RA = [42:45]';
    
    LTIB_CS = quat2rotm(x(idx_q_LA)');
    RTIB_CS = quat2rotm(x(idx_q_RA)');
    PELV_CS = quat2rotm(x(idx_q_MP)');
    
    LKNE = x(idx_pos_LA,1) + d_ltibia*LTIB_CS(:,3);
    RKNE = x(idx_pos_RA,1) + d_rtibia*RTIB_CS(:,3);
    LFEP = x(idx_pos_MP,1) + d_pelvis/2*PELV_CS(:,2);
    RFEP = x(idx_pos_MP,1) - d_pelvis/2*PELV_CS(:,2);
    
    LFEM_z = x(idx_pos_MP,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idx_pos_MP,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;
    
    LFEM_z = LFEM_z/norm(LFEM_z);
    RFEM_z = RFEM_z/norm(RFEM_z);
    
    LFEM_y = LTIB_CS(:,2);
    RFEM_y = RTIB_CS(:,2);
    
    LFEM_x = cross(LFEM_y,LFEM_z);
    RFEM_x = cross(RFEM_y,RFEM_z);
    
    LFEM_CS = [LFEM_x, LFEM_y, LFEM_z];
    RFEM_CS = [RFEM_x, RFEM_y, RFEM_z];
    
    qLFEM = rotm2quat(LFEM_CS)';
    qRFEM = rotm2quat(RFEM_CS)';
    
    x_rec(:,k) = x;
    xa_rec(:,k) = [LFEP; LKNE; RFEP; RKNE];
    qFEM(:,k) = [qLFEM; qRFEM];
end
end

function [x1,X_post,Px,X_dev,S]=ut(g,X,Wm,Wc,nStates,Q,fs)

%Unscented propogation Transformation
%Input:
%        g: nonlinear map
%        X: sigma points
%       Wm: weights for mean
%       Wc: weights for covraiance
%        n: number of outputs of f
%        Q: additive covariance
%       fs: sampleing frequency
%Output:
%        x1: transformed mean
%    X_post: transformed smapling points
%        Px: transformed covariance
%     X_dev: transformed deviations
%         S: "square root" of covariance (col. factor)

L = size(X,2);
x1 = zeros(nStates,1);
X_post = zeros(nStates,L);
for k=1:L
    X_post(:,k)=g(X(:,k),fs);
    x1=x1+Wm(k)*X_post(:,k);
end
X_dev=X_post-x1(:,ones(1,L));
Px=X_dev*diag(Wc)*X_dev'+Q;
A_tmp = [sqrt(Wc(2))*(X_post(:,2:L)-x1) sqrt(Q)]';
if issparse(A_tmp)
    S = qr(A_tmp,0);
else
    S = triu(qr(A_tmp,0));
end
S = S(1:length(x1),1:length(x1));
if Wc(1) < 0
    S = cholupdate(S, (X_post(:,1)-x1), '-');
else
    S = cholupdate(S, (X_post(:,1)-x1), '+');
end
end

function X=genSig(x,S,mu)
%Sigma points symetrically distributed around initial point
%Inputs:
%       x: reference point
%       P: covariance
%       mu: coefficient
%Output:
%       X: Sigma points

% A = mu*chol(P)'; %cholesky factorization
Y = x(:,ones(1,length(x)));
X = [x Y+mu*S' Y-mu*S'];
end
%for state vec with pos vel accel, quat, ang.vel
function [xp] = f_pvaqw(x,fs)
dt = 1/fs;
dt2 = 1/2*(1/fs)^2;
As(1:3,:) = [eye(3,3) dt*eye(3,3) dt2*eye(3,3) zeros(3,7)];
As(4:6,:) = [zeros(3,3) eye(3,3) dt*eye(3,3) zeros(3,7)];
As(7:9,:) = [zeros(3,6) eye(3,3) zeros(3,7)];
As(10:16,:) = [zeros(7,9) eye(7,7)];
Af = [As zeros(16,32); zeros(16,16) As zeros(16,16); zeros(16,32) As];
xp = Af*x;
%mp
e0 = xp(10); e1 = xp(11); e2 = xp(12); e3 = xp(13);
epm = [-e1 -e2 -e3;...
    e0 -e3 e2;...
    e3 e0 -e1;...
    -e2 e1 e0];
xp(10:13) = xp(10:13) + 0.5*epm*xp(14:16)*dt;
%LA
e0 = xp(26); e1 = xp(27); e2 = xp(28); e3 = xp(29);
epm = [-e1 -e2 -e3;...
    e0 -e3 e2;...
    e3 e0 -e1;...
    -e2 e1 e0];
xp(26:29) = xp(26:29) + 0.5*epm*xp(30:32)*dt;
%RA
e0 = xp(42); e1 = xp(43); e2 = xp(44); e3 = xp(45);
epm = [-e1 -e2 -e3;...
    e0 -e3 e2;...
    e3 e0 -e1;...
    -e2 e1 e0];
xp(42:45) = xp(42:45) + 0.5*epm*xp(46:48)*dt;
end

%for state vec with pos vel accel, quat, ang.vel
function [hp] = h_pvaqw(x,fs)
Hs(1:3,:) = [zeros(3,6) eye(3,3) zeros(3,7)];
Hs(4:7,:) = [zeros(4,9) eye(4,4) zeros(4,3)];
Hs(8:10,:) = [zeros(3,13) eye(3,3)];
Hf = [Hs zeros(10,32); zeros(10,16) Hs zeros(10,16); zeros(10,32) Hs];
hp = Hf*x;
end

function y = L2Dist(x,x0,S)
%x^2 is monotomically increasing at any point not at 0, so traditional
%L2 norm involving sqrt is unnecessary, can use x^2 to find same
%location of min cost in constrained region with less computational
%cost.
%using inverse of covariance rather than I will make state est over time
%less smooth but more accurate over the average of the interval

%y = (x-x0)'*(x-x0);
% add index specifying pos. of mp la and ra rather than including q
qW = 1;%norm(x0); %weighting for q deviation to account for diffin units between pos and q.
%need to adapt x to rel pos of LA and RA to prevent large data recording
%problems
%posIdx = [];
res = (x-x0);
qIdx = [10:13 26:29 42:45]';
res(qIdx) = qW*res(qIdx);
res = S\res; %add res*inv(S) to scale cost by certainty
n = 10; %have also tried 2,4,6,8,14,100,1000 around >= 14 greatly increases speed of finding solution, not much difference etween 100 and 1000
%y = (res'*res)^n;
y = sum(res.^n);
end

function [c, ceq] = hingeJoint_constrNL_q(xhat,...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia)
%% Pos-level knee joint constraint
idx_pos_MP = [1:3]';
idx_pos_LA = [17:19]';
idx_pos_RA = [33:35]';

idx_vel_LA = [20:22]';
idx_vel_RA = [36:38]';

idx_q_MP = [10:13]';
idx_q_LA = [26:29]';
idx_q_RA = [42:45]';

I_N = eye(length(xhat));
%       pos     vel         accel       quat    ang.vel      pos     vel         accel       quat    ang.vel        pos     vel         accel       quat    ang.vel
D = [-eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3);
    -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4)  zeros(3,3)];

%rotate from sensor to gfr frame
LTIB_CS = quat2rotm(xhat(idx_q_LA,1)');
RTIB_CS = quat2rotm(xhat(idx_q_RA,1)');
PELV_CS = quat2rotm(xhat(idx_q_MP,1)');


LKNE = xhat(idx_pos_LA,1) + d_ltibia*LTIB_CS(:,3);
RKNE = xhat(idx_pos_RA,1) + d_rtibia*RTIB_CS(:,3);

% calculate the z axis of the femur
LFEM_z = xhat(idx_pos_MP,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
RFEM_z = xhat(idx_pos_MP,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;

% normalize z-axis of femur
LFEM_z = LFEM_z/norm(LFEM_z);
RFEM_z = RFEM_z/norm(RFEM_z);

% calculate the z axis of the tibia
LTIB_z = LTIB_CS(:,3);
RTIB_z = RTIB_CS(:,3);

% calculate alpha_lk and alpha_rk
%alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
%alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));

%calculate alpha_lk and alpha_rk using atan2
LFEM_z__N = LFEM_z;
RFEM_z__N = RFEM_z;

% _TIB_CS is _TIB_CS described in wrod frame, or rotm from tib2world frame
% therefore, inverse is from wrold to tib frame
LFEM_z__TIB = LTIB_CS\LFEM_z__N;
RFEM_z__TIB = RTIB_CS\RFEM_z__N;

alpha_lk = atan2(-LFEM_z__TIB(1),LFEM_z__TIB(3));
alpha_rk = atan2(-RFEM_z__TIB(1),RFEM_z__TIB(3));

% setup the constraint equations
d_k = [ (d_pelvis/2*PELV_CS(:,2) ...
    -d_lfemur*cos(alpha_lk)*LTIB_CS(:,3) ...
    +d_lfemur*sin(alpha_lk)*LTIB_CS(:,1) ...
    -d_ltibia*LTIB_CS(:,3)) ; ...
    (-d_pelvis/2*PELV_CS(:,2)+ ...
    -d_rfemur*cos(alpha_rk)*RTIB_CS(:,3) ...
    +d_rfemur*sin(alpha_rk)*RTIB_CS(:,1) ...
    -d_rtibia*RTIB_CS(:,3)) ];

%% Vel-Level Knee-Joint Constraint
for i = 1:2
    %% Variable assignment based on knee cosntructing constraints
    %(left or right knee)
    
    l_pel = d_pelvis;
    
    wxPEL = xhat(14);
    wyPEL = xhat(15);
    wzPEL = xhat(16);
    
    q1PEL = xhat(10);
    q2PEL = xhat(11);
    q3PEL = xhat(12);
    q4PEL = xhat(13);
    
    dxv = [xhat(20) - xhat(4); %left ankle_pel rel vel in Nx
        xhat(21) - xhat(5); % || in Ny
        xhat(22) - xhat(6);%  || in Nz
        xhat(36) - xhat(4);% right ankle _pel rel vel in Nx
        xhat(37) - xhat(5);%    || Ny
        xhat(38) - xhat(6)];%   || Nz
    
    if i == 1 %left knee
        
        l_tib = d_ltibia;
        l_fem = d_lfemur;
        
        qKN = alpha_lk;
        
        wxTIB = xhat(30);
        wyTIB = xhat(31);
        wzTIB = xhat(32);
        
        q1TIB = xhat(26);
        q2TIB = xhat(27);
        q3TIB = xhat(28);
        q4TIB = xhat(29);
        
    elseif i == 2 %right knee
        l_tib = d_rtibia;
        l_fem = d_rfemur;
        
        qKN = alpha_rk;
        
        wxTIB = xhat(46);
        wyTIB = xhat(47);
        wzTIB = xhat(48);
        
        q1TIB = xhat(42);
        q2TIB = xhat(43);
        q3TIB = xhat(44);
        q4TIB = xhat(45);
    else
        error('Knee not properly specifid for vel. cosntraint')
    end
    
    %% constraint calc
    %qKN = -qKN;
    %w_KN_scl is a scalar value of the knee joint angular velocity
    w_KN_scl = wyTIB + 2*(q1TIB*q4TIB-q2TIB*q3TIB)*(wxTIB*cos(qKN)+...
        wzTIB*sin(qKN))/(2*sin(qKN)*(q1TIB*q3TIB+q2TIB*q4TIB)+...
        cos(qKN)*(-1+2*q1TIB^2+2*q2TIB^2));
    %w_KN_scl = -w_KN_scl;
    %relVel_ANK_PELo_N is the relative velocities of the ankle from the pelvis,
    %described in frame N (eqivalent to gfr)
    relVel_ANK_PELo_N = [(l_pel*wxPEL*(q1PEL*q3PEL+q2PEL*q4PEL)-...
        l_tib*wyTIB*(-1+2*q1TIB^2+2*q2TIB^2)-0.5*l_pel*wzPEL*(-1+2*q1PEL^2+2*q2PEL^2)...
        -2*(q1TIB*q4TIB-q2TIB*q3TIB)*(l_tib*wxTIB+l_fem*(wxTIB*cos(qKN)+wzTIB*sin(qKN)))...
        -l_fem*(2*sin(qKN)*(q1TIB*q3TIB+q2TIB*q4TIB)+cos(qKN)*(-1+2*q1TIB^2+...
        2*q2TIB^2))*(wyTIB-w_KN_scl));...
        ((-1+2*q1TIB^2+2*q3TIB^2)*(l_tib*wxTIB+l_fem*(wxTIB*cos(qKN)+wzTIB*sin(qKN)))...
        -2*l_tib*wyTIB*(q1TIB*q4TIB+q2TIB*q3TIB)-l_pel*wzPEL*(q1PEL*q4PEL+q2PEL*q3PEL)...
        -l_pel*wxPEL*(q1PEL*q2PEL-q3PEL*q4PEL)-2*l_fem*(cos(qKN)*(q1TIB*q4TIB+q2TIB*q3TIB)...
        -sin(qKN)*(q1TIB*q2TIB-q3TIB*q4TIB))*(wyTIB-w_KN_scl));...
        (l_pel*wzPEL*(q1PEL*q3PEL-q2PEL*q4PEL)+2*l_tib*wyTIB*(q1TIB*q3TIB-q2TIB*q4TIB)...
        +0.5*l_pel*wxPEL*(-1+2*q1PEL^2+2*q4PEL^2)+2*(q1TIB*q2TIB+...
        q3TIB*q4TIB)*(l_tib*wxTIB+l_fem*(wxTIB*cos(qKN)+wzTIB*sin(qKN)))+...
        l_fem*(2*cos(qKN)*(q1TIB*q3TIB-q2TIB*q4TIB)-sin(qKN)*(-1+2*q1TIB^2 ...
        +2*q4TIB^2))*(wyTIB-w_KN_scl))];
    
    %left, then right knee constrained
    if i == 1
        d_k_v(1:3,1) = relVel_ANK_PELo_N;
    elseif i == 2
        d_k_v(4:6,1) = relVel_ANK_PELo_N;
    end
end
%% construct constraint for fmincon format
% res = d_k - D*xhat;
% ceq = res'*res;
kv = 0.1;
pRes = (d_k - D*xhat);
vRes = (d_k_v - dxv);
kq = 100;
qMPRes = (xhat(idx_q_MP)'*xhat(idx_q_MP))-1;
qLARes = (xhat(idx_q_LA)'*xhat(idx_q_LA))-1;
qRARes = (xhat(idx_q_RA)'*xhat(idx_q_RA))-1;

ceq = [ pRes + kv*vRes;
        kq*qMPRes;
        kq*qLARes;
        kq*qRARes];
    

maxFootVel = 12.5; %m/s

% c = [norm(xhat(idx_vel_LA))-maxFootVel;
%     norm(xhat(idx_vel_RA))-maxFootVel];
qSafetyFactor = deg2rad(3);
c = [-alpha_lk;
    (alpha_lk+qSafetyFactor) - 0.5*pi;
    -alpha_rk;
    (alpha_rk+qSafetyFactor) - 0.5*pi];
%TODO (or at least consider):
% add quaternion constraints (motion, vel. accel) as in 331 text
% add accel-level knee-joint constraint
% Tune weighting on pos-vel knee constraints (kd)
% add additional physiological params (ang vel, etc.)
% Add ZVUPT
% Add trailing accel rec, to limit jerk profile->smooth/limit accel
end
