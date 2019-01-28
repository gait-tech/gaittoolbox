%% Explore circle constraint toy example

x = [1.5*cosd(45) 1.5*sind(45)]';
P = 1e-2*[1 0;
         0 1];
dLFemur = 1.0;
nSamples = 500;

Pbuf = (P+P.')/2;
r = mvnrnd(x, Pbuf, nSamples);

r2 = zeros(nSamples, 2);
for i=1:nSamples
    [newx, newP] = applySCKF(r(i,:)', P, dLFemur);
    r2(i, :) = newx;
end

% manual sckf
sckfAlpha = 0.1;
sckfThreshold = 100;
sckfMaxIter = 50;
applyCstr = 1;
I_N = eye(2, 2);

x_tilde = x;
P_tilde = P;

for i=0:sckfMaxIter
    % calculate the z axis of femur and tibia
    LFEM_z = x_tilde(:, 1);
    g_dlfem = norm(LFEM_z, 2);
    LFEM_z2 = x_tilde(:, 1) - [0; 1];
    g_dlfem2 = norm(LFEM_z2, 2);


    D = zeros(2, 2);
    D(1, 1:2) = LFEM_z'/g_dlfem;
    D(2, 1:2) = LFEM_z2'/g_dlfem2;
    res = [dLFemur - g_dlfem; 
           dLFemur - g_dlfem2];
    if i==0
        R0 = sckfAlpha*D*P_tilde*D';
    end

    Ri = R0*exp(-i);
    Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ diag(D * P_tilde* D');
    Si(isnan(Si)) = sckfThreshold+1;
    if sum(Si < sckfThreshold) == 0, break, end

    switch mod(applyCstr, 10)
        case 1
            Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
            P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
        case 2
            Kk = D'*(D*D'+Ri)^(-1);
            P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
        case 3
            Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
        case 4
            Kk = D'*(D*D'+Ri)^(-1);
        case 5
            Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
            P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
        case 6
            Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
        case 7
            Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
            P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
        case 8
            Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
        otherwise
            Kk = 0;
    end

    d = res + D*x_tilde;
    dx = Kk*(res);
    x_tilde = x_tilde + dx;

    if true
        clf; hold on; axis square;
        viscircles([0 0], 1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.1);
        viscircles([0 1], 1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.1);

        scatter(r(:, 1), r(:, 2), 'b.');
        scatter(r2(:, 1), r2(:, 2), 'r.');
        line([d(1) / D(1, 1) 0], [0 d(1) / D(1, 2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
        line([(d(2)-D(2,2)*2) / D(2, 1) 0], [2 d(2) / D(2, 2)], 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1);
        Pbuf = (P_tilde+P_tilde.')/2;
        r3 = mvnrnd(x_tilde, Pbuf, nSamples);
        scatter(r3(:, 1), r3(:, 2), 'g.');
        legend({'unconstrained', 'cstr est cov sample', 'cstr act cov sample'});
        xlim([0 2]); ylim([0 2]);
        pause(1);
    end
end

cstrx = x_tilde;
cstrP = P_tilde;
        
function [newx, newP, D, d] = applySCKF(x, P, dLFemur)
    sckfAlpha = 0.1;
    sckfThreshold = 100;
    sckfMaxIter = 50;
    applyCstr = 1;
    I_N = eye(2, 2);
    
    x_tilde = x;
    P_tilde = P;
            
    for i=0:sckfMaxIter
        % calculate the z axis of femur and tibia
        LFEM_z = x_tilde(:, 1);
        g_dlfem = norm(LFEM_z, 2);
        LFEM_z2 = x_tilde(:, 1) - [0; 1];
        g_dlfem2 = norm(LFEM_z2, 2);

        
        D = zeros(2, 2);
        D(1, 1:2) = LFEM_z'/g_dlfem;
        D(2, 1:2) = LFEM_z2'/g_dlfem2;
        res = [dLFemur - g_dlfem; 
               dLFemur - g_dlfem2];
        if i==0
            R0 = sckfAlpha*D*P_tilde*D';
        end

        Ri = R0*exp(-i);
        Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ diag(D * P_tilde* D');
        Si(isnan(Si)) = sckfThreshold+1;
        if sum(Si < sckfThreshold) == 0, break, end

        switch mod(applyCstr, 10)
            case 1
                Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
            case 2
                Kk = D'*(D*D'+Ri)^(-1);
                P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
            case 3
                Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
            case 4
                Kk = D'*(D*D'+Ri)^(-1);
            case 5
                Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
            case 6
                Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
            case 7
                Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
            case 8
                Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
            otherwise
                Kk = 0;
        end

        d = res + D*x_tilde;
        dx = Kk*(res);
        x_tilde = x_tilde + dx;
    end
    
    newx = x_tilde;
    newP = P_tilde;
end