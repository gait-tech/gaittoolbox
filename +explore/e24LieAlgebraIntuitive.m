addpath('liese3lib');

C0 = eul2rotm([3/4*pi, pi/2, pi]);
r0 = [0 0 0]';
T0 = [C0 r0; 0 0 0 1];

C1 = eul2rotm([1/4*pi, 1/4*pi, 1/3*pi]);
w0 = rot2vec(C0'*C1)';

clf;
pelib.viz.plotPosOri(T0(1:3,4), T0(1:3,1:3));

% for i=w
    xi = [1 1 1 w0]';

    B_R_W = T0(1:3,1:3)';
    xiOld = xi;

    xi(4:6) = xi(4:6);
    J = vec2jac(xi(4:6));
    xi(1:3) = inv(J)*B_R_W*xi(1:3);
    bigxi = vec2tran(xi);

    res0 = T0*bigxi;
    res1 = T0(1:3,4) + xiOld(1:3);
    
    fprintf('Lie: %d %d %d \nVec: %d %d %d\n', res0(1:3,4)', res1')
    pelib.viz.plotPosOri(res0(1:3,4), res0(1:3,1:3));
    scatter3(res1(1), res1(2), res1(3));
% end
xlabel('x'); ylabel('y'); zlabel('z'); axis square;