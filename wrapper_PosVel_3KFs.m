function [ xHat,Px,yHat,Py,zHat,Pz ] = wrapper_PosVel_3KFs( fs,acc,x0,sigma_acc,bUpdate)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% x0 is [3 x 1] containing the initial position

xKF = KF_PosVel('x0',[x0(1) x0(4)]','fs',fs,'sigma_acc',sigma_acc);
yKF = KF_PosVel('x0',[x0(2) x0(5)]','fs',fs,'sigma_acc',sigma_acc);
zKF = KF_PosVel('x0',[x0(3) x0(6)]','fs',fs,'sigma_acc',sigma_acc);

[N_ACC,~] = size(acc);
xHat = nan(N_ACC,2);
yHat = nan(N_ACC,2);
zHat = nan(N_ACC,2);

Px = nan(2,2,N_ACC);
Py = nan(2,2,N_ACC);
Pz = nan(2,2,N_ACC);

sigma_vel = 0.1;
for n = 1:N_ACC
% ---- Prediction Steps ----    
    % x position
    xKF.PredictStep(acc(n,1));
    xHat(n,:) = xKF.xHat;
    Px(:,:,n) = xKF.P;
    % y position
    yKF.PredictStep(acc(n,2));
    yHat(n,:) = yKF.xHat;
    Py(:,:,n) = yKF.P;
    % z position
    zKF.PredictStep(acc(n,3));
    zHat(n,:) = zKF.xHat;
    Pz(:,:,n) = zKF.P;
% ---- Update Steps ----    
    if bUpdate(n,1)
        xKF.UpdateStep(sigma_vel);
        xHat(n,:) = xKF.xHat;
        Px(:,:,n) = xKF.P;
    end
    if bUpdate(n,2)
        yKF.UpdateStep(sigma_vel);
        yHat(n,:) = yKF.xHat;
        Py(:,:,n) = yKF.P;
    end
    if bUpdate(n,3)
        zKF.UpdateStep(sigma_vel);
        zHat(n,:) = zKF.xHat;
        Pz(:,:,n) = zKF.P;
    end    
end

end

