%cukf_v6 (included quaternions in state, first using last timestep as 
%prediction for next, v7 will change this implmentation and remove accel 
%from state ve) still constaints implemented knee pinjoint constraint
%on SR-UKF with nonlin projection of estimated x onto constrained space
function [ x_rec, xa_rec, qFEM ] = cukf_v6(x,P,Q,R,N_MP,nMeas,acc,fs,...
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
xa_rec = nan(12,N_MP); % 12 = femur proximal and distal points, both legs
qFEM = nan(8,N_MP); %8 = 4 per quaterinoin and 2 legs
%vanderwerve 2001 paper wm, wc, and constant definitions

for k=1:N_MP
    disp('k')
    disp(k)
    z = [acc(k,1:3)'; q_MP(k,:)'; acc(k,4:6)'; q_LA(k,:)'; acc(k,7:9)'; q_RA(k,:)'];
    %could constrain raw measurements here
    
    X = genSig(x,S,mu);                            %sigma points around x
    
    [x1,X_post,Px,X_dev,S] = ut(@f_waaq,X,Wm,Wc,L,Q,fs);          %unscented transformation of process
    
    [z1,Z_post,Pz,Z_dev,Sy] = ut(@h_waaq,X_post,Wm,Wc,m,R,fs);       %unscented transformation of measurments
    %constrain measuremnts
    Pxz=X_dev*diag(Wc)*Z_dev';                        %transformed cross-covariance
    
    K=(Pxz/Sy)/Sy';
    %K=Pxz/(Pz);
%     disp('dim K')
%     disp(size(K));
    P = Px-K*Pxz';
if isConstr
        xhat = x1;
        options = optimoptions('fmincon','Algorithm','interior-point','Display','off',...
            'OptimalityTolerance', 1e-4, 'ConstraintTolerance', 1e-4,...
            'MaxFunctionEvaluations',15000); % run interior-point algorithm
    x1 = fmincon(@(x1) L2Dist(x1,xhat,S),xhat,[],[],[],[],[],[],@(x1) hingeJoint_constrNL_q(x1,...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia),options);
end
    x=x1+K*(z-z1);                              %measurement update
    
%     if isConstr
%         xhat = x;
%         options = optimoptions('fmincon','Algorithm','sqp','Display','off',...
%             'OptimalityTolerance', 1e-4, 'ConstraintTolerance', 1e-4,...
%             'MaxFunctionEvaluations',15000); % run interior-point algorithm
%     x = fmincon(@(x) L2Dist(x,xhat),xhat,[],[],[],[],[],[],@(x) hingeJoint_constrNL_q(x,...
%     d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia),options);

%         [x] = hingeJoint_constrNL(x,hjc,P,k,q_MP, q_LA, q_RA,...
%             d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia);
%     end
    U = K*Sy';
    for i = 1:m
        S = cholupdate(S, U(:,i), '-');
    end
    %% gen augmented x return values
    idx_pos_MP = [1:3]';
    idx_pos_LA = [14:16]';
    idx_pos_RA = [27:29]';
    idx_q_MP = [10:13]';
    idx_q_LA = [23:26]';
    idx_q_RA = [36:39]';
    
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

function [xp] = f_waaq(x,fs)
dt = 1/fs;
dt2 = 1/2*(1/fs)^2;
As(1:3,:) = [eye(3,3) dt*eye(3,3) dt2*eye(3,3) zeros(3,4)];
As(4:6,:) = [zeros(3,3) eye(3,3) dt*eye(3,3) zeros(3,4)];
As(7:13,:) = [zeros(7,6) eye(7,7)];
Af = [As zeros(13,26); zeros(13,13) As zeros(13,13); zeros(13,26) As];
xp = Af*x;
end

function [hp] = h_waaq(x,fs)
Hs(1:3,:) = [zeros(3,6) eye(3,3) zeros(3,4)];
Hs(4:7,:) = [zeros(4,9) eye(4,4)];
Hf = [Hs zeros(7,26); zeros(7,13) Hs zeros(7,13); zeros(7,26) Hs];
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
qW = norm(x0); %weighting for q deviation to account for diffin units between pos and q.
%need to adapt x to rel pos of LA and RA to prevent large data recording
%problems
%posIdx = [];
res = (x-x0);
qIdx = [10:13 23:26 36:39]';
res(qIdx) = qW*res(qIdx);
res = S\res; %add res*inv(S) to scale cost by certainty
n = 1000; %have also tried 2,4,6,8,14,100,1000 around >= 14 greatly increases speed of finding solution, not much difference etween 100 and 1000
y = (res'*res)^n;

end

function [c, ceq] = hingeJoint_constrNL_q(xhat,...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia)
%   if hingejoint_constraint > 0 && hingejoint_constraint <= 8
% calculate the location of the knee
idx_pos_MP = [1:3]';
idx_pos_LA = [14:16]';
idx_pos_RA = [27:29]';
idx_q_MP = [10:13]';
idx_q_LA = [23:26]';
idx_q_RA = [36:39]';

I_N = eye(length(xhat));

D = [-eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,4);
    -eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4) zeros(3,3) zeros(3,3) zeros(3,3) zeros(3,4) eye(3,3) zeros(3,3) zeros(3,3) zeros(3,4)];

LTIB_CS = quat2rotm(xhat(idx_q_LA,1)');
RTIB_CS = quat2rotm(xhat(idx_q_RA,1)');
PELV_CS = quat2rotm(xhat(idx_q_MP,1)');

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

% res = d_k - D*xhat;
% ceq = res'*res;

ceq = d_k - D*xhat;

c = [];
%note: still want to use kalman gain to scale res or costfcn?
%weight cost function by inverse covariance matrix
%after 3q. 28 in simon 2010
%dx = Kk*(res);

end
