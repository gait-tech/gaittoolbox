%% Explore femur constraint toy example

x = [-0.5*cosd(45) 0.5*sind(45) 0 0]';
P = 1e-2*[1 0 0 0;
         0 1 0 0;
         0 0 1 0;
         0 0 0 1];
H = [0 0 1 0;
     0 0 0 1];
R = [1 0; 0 1e-5];

for i=1:500
%     draw(x, P, 500);
%     pause(1/1000);
    
    [newx, newP] = applySCKF(x, P, 1.0);
    x = newx;
%     P = newP;
    norm(x(1:2)-x(3:4))
    
    draw(x, P, 500); title('Thigh Constraint Part B');
    pause(1/1000);
    % fake prediction update that reduces length
    v = x(1:2) - x(3:4);
    v = v / norm(v);
    x = x - 0.1*[v; 0; 0];
    P = P + 1e-5*eye(4,4);
    
    % measurement update
    res = [0; 0] - H*x;
    K = P * H' /(H*P*H' + R);
    x = x + K * res;
    P = (eye(4,4) - K * H) * P;
end

function [newx, newP] = applySCKF(x, P, dLFemur)
    sckfAlpha = 0.1;
    sckfThreshold = 100;
    sckfMaxIter = 50;
    applyCstr = 1;
    I_N = eye(4, 4);
    
    x_tilde = x;
    P_tilde = P;
            
    for i=0:sckfMaxIter
        % calculate the z axis of femur and tibia
        LFEM_z = x_tilde(1:2, 1) - x_tilde(3:4, 1);
        g_dlfem = norm(LFEM_z, 2);

        D = zeros(1, 4);
        D(1, 1:2) = LFEM_z'/g_dlfem;
        D(1, 3:4) = -LFEM_z'/g_dlfem;
        res = [dLFemur - g_dlfem];
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

        dx = Kk*(res);
        x_tilde = x_tilde + dx;
    end
    
    newx = x_tilde;
    newP = P_tilde;
end

function draw(x, P, nSamples)
    clf; hold on; axis square;
    line(x([1 3]), x([2 4]), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1); 
    line([x(3) x(3)], [x(4) x(4)-1], 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1); 
    viscircles(x(3:4,1)', 1, 'Color', 'k', 'LineStyle', '--', 'LineWidth', 0.1);
    
    Pbuf = P(1:2,1:2);
    Pbuf = (Pbuf+Pbuf.')/2;
    r = mvnrnd(x(1:2)', Pbuf, nSamples);
    scatter(r(:, 1), r(:, 2), 'bx');
    
    Pbuf = P(3:4,3:4);
    Pbuf = (Pbuf+Pbuf.')/2;
    r = mvnrnd(x(3:4)', Pbuf, nSamples);
    scatter(r(:, 1), r(:, 2), 'gx');
    
    prob = mvnpdf(r, x(3:4)', Pbuf);
    prob = repelem(prob, 180, 1);
    r2 = zeros(nSamples*180, 2);
    for i=1:nSamples
        for j=1:180
            r2((i-1)*180+j, :) = r(i, :) + [-cosd(j-90) sind(j-90)];
        end
    end
    idx = randsample(1:(nSamples*180), nSamples, true, prob);
    r = r2(idx, :);
    scatter(r(:, 1), r(:, 2), 'b.');

    r2 = mean(r);
    scatter(r2(:, 1), r2(:, 2), 'bo');
%     r = mvnrnd(x, P, nSamples);
%     r2 = zeros(nSamples, 4);
%     P2 = eye(4,4); P2(3:4,3:4) = 0;
%     for i=1:nSamples
%         r2(i,:) = applySCKF(r(i,:)', P2, 1);
%     end
%     scatter(r2(:, 1), r2(:, 2), 'y^');
%     scatter(r2(:, 3), r2(:, 4), 'k^');
    
    theta = mod(270-atan2d(x(2)-x(4), x(1)-x(3)), 360);
    text(x(3)-0.5, x(4), sprintf('%.2f deg', theta));
    scatter(x([1 3]), x([2 4]), 'k*');
    legend({'Thigh', 'Shanks', ...
        'Pelv Pos Cov Est (lin)', 'Knee Pos Cov Est', ...
        'Pelv Pos Cov Act (nonlin)', 'Pelvis Pos Mean Act', ...
        'Actual Pelv and Knee Pos'});
end