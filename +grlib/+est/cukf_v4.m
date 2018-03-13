%cukf_v4 implemented knee pinjoint constraint on SR-UKF
function [ x_rec, xa_rec, qFEM ] = cukf_v4(x,P,Q,R,N_MP,nMeas,acc,fs,...
    q_MP, q_LA, q_RA, d_pelvis, d_lfemur, d_rfemur,...
    d_ltibia, d_rtibia,isConstr)
% Unscented Kalman filter for state estimation
% UKF   Unscented Kalman Filter for nonlinear dynamic systems
% [x, P] = ukf(f,x,P,h,z,Q,R) returns state estimate, x and state covariance, P
% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
%           x_k+1 = f(x_k) + w_k
%           z_k   = h(x_k) + v_k
% where w ~ N(0,Q) meaning w is gaussian noise with covariance Q
%       v ~ N(0,R) meaning v is gaussian noise with covariance R
% Inputs:   f: function handle for f(x)
%           x: "a priori" state estimate
%           P: "a priori" estimated state covariance
%           h: fanction handle for h(x)
%           z: current measurement
%           Q: process noise covariance
%           R: measurement noise covariance
% Output:   x: "a posteriori" state estimate
%           P: "a posteriori" state covariance

L = length(x);                                 %numer of states
m = nMeas;
alpha = 1e-4;     %should be small (0<alpha<1)                            %default, tunable
ki = 0;            %start with 0 change if needed                           %default, tunable
beta = 2;     % assuming gaussian prior                                %default, tunable
hjc = 4;%or 3 - see hingeJoint_cosntr func. for details
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
xa_rec = nan(12,N_MP);
qFEM = nan(8,N_MP);
%vanderwerve 2001 paper wm, wc, and constant definitions

for k=1:N_MP
    
    z = acc(k,:)';
    %could constrain raw measurements here
    
    X = genSig(x,S,mu);                            %sigma points around x
    
    [x1,X_post,Px,X_dev,S] = ut(@f,X,Wm,Wc,L,Q,fs);          %unscented transformation of process
    
    [z1,Z_post,Pz,Z_dev,Sy] = ut(@h,X_post,Wm,Wc,m,R,fs);       %unscented transformation of measurments
    %constrain measuremnts
    Pxz=X_dev*diag(Wc)*Z_dev';                        %transformed cross-covariance
    
    K=(Pxz/Sy)/Sy';
    %K=Pxz/(Pz);
    P = Px-K*Pxz';
    if isConstr
        [x1] = hingeJoint_constr(x1,hjc,P,k,q_MP, q_LA, q_RA,...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia);
    end
    x=x1+K*(z-z1);                              %measurement update
    if isConstr
        [x] = hingeJoint_constr(x,hjc,P,k,q_MP, q_LA, q_RA,...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia);
    end
    U = K*Sy';
    for i = 1:m
        S = cholupdate(S, U(:,i), '-');
    end
    
    idx_pos_MP = [1:3]';
    idx_pos_LA = [10:12]';
    idx_pos_RA = [19:21]';
    
    LTIB_CS = quat2rotm(q_LA(k,:));
    RTIB_CS = quat2rotm(q_RA(k,:));
    PELV_CS = quat2rotm(q_MP(k,:));
    
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
%        R: additive covariance
%Output:
%        y: transformed mean
%        Y: transformed smapling points
%        P: transformed covariance
%       Y1: transformed deviations

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

function [x_constr] = constr(x,P,n,q_MP, q_LA, q_RA,...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia)
% Coplanar pos constraint (9 el. state per sensor)
D = [-eye(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
    -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3)];

%roation matrices from body frame to ground frame
LTIB_CS = quat2rotm(q_LA(n,:));
RTIB_CS = quat2rotm(q_RA(n,:));
PELV_CS = quat2rotm(q_MP(n,:));

% calculate the location of the knee
LKNE = x(10:12,1) + d_ltibia*LTIB_CS(:,3);
RKNE = x(19:21,1) + d_rtibia*RTIB_CS(:,3);

% calculate the z axis of the femur
LFEM_z = x(1:3,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
RFEM_z = x(1:3,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;

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

% res = d_k-D*x;
% Kk = P*D'/(D*P*D');
% dx = Kk*(res);
dx = P*D'/(D*P*D')*(d_k-D*x);
x_constr = x + dx;
end

function [x_constr] = hingeJoint_constr(xhat,hjc,P,n,q_MP, q_LA, q_RA,...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia)
%   if hingejoint_constraint > 0 && hingejoint_constraint <= 8
% calculate the location of the knee
idx_pos_MP = [1:3]';
idx_pos_LA = [10:12]';
idx_pos_RA = [19:21]';

I_N = eye(length(xhat));

D = [-eye(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3);
    -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,3) eye(3,3) zeros(3,3) zeros(3,3)];

LTIB_CS = quat2rotm(q_LA(n,:));
RTIB_CS = quat2rotm(q_RA(n,:));
PELV_CS = quat2rotm(q_MP(n,:));

LKNE = xhat(idx_pos_LA,1) + d_ltibia*LTIB_CS(:,3);
RKNE = xhat(idx_pos_RA,1) + d_rtibia*RTIB_CS(:,3);

% calculate the z axis of the femur
LFEM_z = xhat(idx_pos_MP,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
RFEM_z = xhat(idx_pos_MP,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;

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

switch (hjc)
%     case 1 % maximum probability estimate + no P update
%         Kk = P*D'*(D*P*D')^(-1);
%     case 2 % least squares estimate + no P update
%         Kk = D'*(D*D')^(-1);
    case 3 % maximum probability estimate + constrained projection + no P update
        Kk = P*D'*(D*P*D')^(-1);
        %A = (I_N-Kk*D)*A;
    case 4 % least squares estimate + constrained projection + no P update
        Kk = D'*(D*D')^(-1);
        %A = (I_N-Kk*D)*A;
%     case 5 % maximum probability estimate
%         Kk = P*D'*(D*P*D')^(-1);
%         P = (I_N-Kk*D)*P;
%     case 6 % least squares estimate
%         Kk = D'*(D*D')^(-1);
%         P = (I_N-Kk*D)*P;
%     case 7 % maximum probability estimate + constrained projection
%         Kk = P*D'*(D*P*D')^(-1);
%         A = (I_N-Kk*D)*A;
%         P = (I_N-Kk*D)*P;
%     case 8 % least squares estimate + constrained projection
%         Kk = D'*(D*D')^(-1);
%         A = (I_N-Kk*D)*A;
%         P = (I_N-Kk*D)*P;
    otherwise
        Kk = 0;
end

dx = Kk*(res);
x_constr = xhat + dx;
end

function [xp] = f(x,fs)
nDim = 3;
nSense = 3;
nStpSns = 9;
nSt = length(x);
xp = zeros(nSt,1);
for i = 0:nStpSns:nSt-nDim
    for j = 1:nSense
        xp(i+j) = x(i+j)+(x(i+j+nDim)*1/fs)+(0.5*(1/fs)^2*x(i+j+2*nDim));...
        xp(i+j+nDim) = x(i+j+nDim)+1/fs*x(i+j+2*nDim);...
        xp(i+j+2*nDim) = x(i+j+2*nDim);
%     disp('i = ')
%     disp(i)
%     disp('j = ')
%     disp(j)
    end
end
%add noise
%xp = xp + diag(Q);
end

function [hp] = h(x,fs)
nDim = 3;
nSense = 3;
nStpSns = 9;
nSt = length(x);
hp = zeros(nDim*nSense,1);
k = 1;
for i = 0:nStpSns:nSt-nDim
    for j = 1:nSense
        hp(k) = x(i+j+2*nDim);
        k = k+1;
    end
end
%add noise
%hp = hp + diag(R);
end