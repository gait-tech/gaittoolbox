addpath('liese3lib');

%% effect of omega on position
% C0 = eul2rotm([3/4*pi, pi/2, pi]);
% r0 = [0 0 0]';
% T0 = [C0 r0; 0 0 0 1];
% 
% C1 = eul2rotm([1/4*pi, 1/4*pi, 1/3*pi]);
% w0 = rot2vec(C0'*C1)';
% 
% clf;
% pelib.viz.plotPosOri(T0(1:3,4), T0(1:3,1:3));
% 
% % for i=w
%     xi = [1 1 1 w0]';
% 
%     B_R_W = T0(1:3,1:3)';
%     xiOld = xi;
% 
%     xi(4:6) = xi(4:6);
%     J = vec2jac(xi(4:6));
%     xi(1:3) = inv(J)*B_R_W*xi(1:3);
%     bigxi = vec2tran(xi);
% 
%     res0 = T0*bigxi;
%     res1 = T0(1:3,4) + xiOld(1:3);
%     
%     fprintf('Lie: %d %d %d \nVec: %d %d %d\n', res0(1:3,4)', res1')
%     pelib.viz.plotPosOri(res0(1:3,4), res0(1:3,1:3));
%     scatter3(res1(1), res1(2), res1(3));
% % end
% xlabel('x'); ylabel('y'); zlabel('z'); axis square;

%% covariance sampling
P0 = [0.1*eye(3,3) zeros(3,3); zeros(3,3) deg2rad(0.1)*eye(3,3)];
P = (P0+P0.')/2;
mu = zeros(6,1);
R0 = eul2rotm(deg2rad([90 0 0]));
p0 = [0 0 1 1]';

nspace = [10 10 10 3 3 3];
% rlim = [ 0.05  0.05  0.05  deg2rad(20)  deg2rad(20)  deg2rad(20); ...
%         -0.05 -0.05 -0.05 -deg2rad(20) -deg2rad(20) -deg2rad(20)];
psigma = 0.1;
osigma = 5;
rlim = [ psigma  psigma  psigma  deg2rad(osigma)  deg2rad(osigma)  deg2rad(osigma); ...
        -psigma -psigma -psigma -deg2rad(osigma) -deg2rad(osigma) -deg2rad(osigma)];
    
[r0, r1, r2, r3, r4, r5] = ndgrid(linspace(rlim(2,1), rlim(1,1), nspace(1)), ...
    linspace(rlim(2,2), rlim(1,2), nspace(2)), ...
    linspace(rlim(2,3), rlim(1,3), nspace(3)), ...
    linspace(rlim(2,4), rlim(1,4), nspace(4)), ...
    linspace(rlim(2,5), rlim(1,5), nspace(5)), ...
    linspace(rlim(2,6), rlim(1,6), nspace(6)) );
r = [r0(:) r1(:) r2(:) r3(:) r4(:) r5(:)];
[nSamples, ~] = size(r);

% SO(3) x R^3
T1s = repelem(eye(4,4),1,1,nSamples);
T1s(1:3,4,:) = r(:,1:3)';
for i=1:nSamples
    T1s(1:3,1:3,i) = R0*vec2rot(r(i,4:6)');
end

p1s = zeros(4, nSamples);
for i=1:nSamples
    p1s(:,i) = T1s(:,:,i)*p0;
end

% SE(3)
T2s = zeros(4,4,nSamples);
T0 = [R0 mu(1:3); zeros(1,3) 1];

for i=1:nSamples
    T2s(:,:,i) = T0*vec2tran(r(i,:)');
end

p2s = zeros(4, nSamples);
for i=1:nSamples
    p2s(:,i) = T2s(:,:,i)*p0;
end

p0s = zeros(3, nSamples);
for i=1:nSamples
    p0s(:,i) = T2s(1:3,1:3,i)*p0(1:3);
end
p0s(4,:) = 1;

% p1s = p1s(:,1);
% p2s = p2s(:,1);
clf; hold on;
sphere;
scatter3(p1s(1,:), p1s(2,:), p1s(3,:), 'r.');
scatter3(p2s(1,:), p2s(2,:), p2s(3,:), 'b.');
scatter3(p0s(1,:), p0s(2,:), p0s(3,:), 'g.');

% pds = p2s-p1s;
% q = quiver3(p1s(1,:),p1s(2,:),p1s(3,:),pds(1,:),pds(2,:),pds(3,:),...
%             '-', 'AutoScale', 'off', 'MaxHeadSize', 0.05, ...
%             'LineWidth', 0.1);
% pds = p1s-p0s;
% q = quiver3(p0s(1,:),p0s(2,:),p0s(3,:),pds(1,:),pds(2,:),pds(3,:),...
%             '--', 'AutoScale', 'off', 'MaxHeadSize', 0.05, ...
%             'LineWidth', 0.5);
axis square;

%% covariance symbolic plot
% syms s t
% slim = [-1 1]; tlim = [-1 1];
% wsigma = deg2rad(90)*ones(3,1);
% vsigma = 1*ones(3,1);
% 
% p1 = wsigma*s + vsigma*t;
% fmesh(p1(1), p1(2), p1(3), [slim tlim], 'Linewidth', 2)
% axis equal
% % 
% % tmin = -1; tmax = 1;
% % 
% % 
% % R0 = eul2rotm(deg2rad([90 0 0]));
% % p0 = [0 1 0 1]';
% % 
% % % SO(3) x R^3
% % p1 = R0*symvec2rot(wsigma*t)*p0(1:3) + vsigma*t;
% % fplot3(p1(1),p1(2),p1(3),[tmin tmax]);
% 
% function [ C ] = symvec2rot( phi )
%     angle = norm(phi);
%     axis = phi/angle;
%     cp = cos(angle);
%     sp = sin(angle);
%     C = cp * eye(3) + (1 - cp) * axis * axis' + sp * symhat(axis);
% end
% 
% function [ vechat ] = symhat( vec )
%     if size(vec,1) == 3  
%         vechat = [  0,     -vec(3),  vec(2);
%                 vec(3),   0    , -vec(1);
%                -vec(2),  vec(1),   0    ];  
%     elseif size(vec,1) == 6
%         vechat = [ symhat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];      
%     end
% end

