dt = 1;
I_N = eye(4);

AxisX = -3:0.1:3;
AxisY = -3:0.1:3;
[GridX, GridY] = meshgrid(AxisX, AxisY);
tIdx1 = 1:2;
tIdx2 = 3:4;

F = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
G = [0.5*dt^2 0; 0 0.5*dt^2; dt 0; 0 dt];
H = [1 1 0 0];

Q = [0.01 0 0 0; 0 0.01 0 0; 0 0 0.01 0; 0 0 0 0.01];
R = [1e-10];

x_pos = [0; 0; 0; 0];
u = [1; 1];
y = [0.5];

P_pos = [0.1 0 0 0; 0 0.2 0 0; 0 0 0.2 0; 0 0 0 0.3];

for i=1:1
    subplot(3, 2, 1); 
    Est_pos1 = mvnpdf([GridX(:) GridY(:)], x_pos(tIdx1)', P_pos(tIdx1,tIdx1));
    Est_pos1 = reshape(Est_pos1, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pos1); grid; 
    title('Aposteriori Position Estimate at t=k');
    xlabel('x'); ylabel('y');

    subplot(3, 2, 2);
    Est_pos2 = mvnpdf([GridX(:) GridY(:)], x_pos(tIdx2)', P_pos(tIdx2,tIdx2));
    Est_pos2 = reshape(Est_pos2, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pos2); grid; title('Aposteriori Velocity Estimate at t=k');
    xlabel('x'); ylabel('y');

    x_pri = F*x_pos+G*u;
    P_pri = F * P_pos * F' + Q;
    
    subplot(3, 2, 3); grid;
    Est_pri1 = mvnpdf([GridX(:) GridY(:)], x_pri(tIdx1)', P_pri(tIdx1,tIdx1));
    Est_pri1 = reshape(Est_pri1, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pri1); grid; title('Apriori Position Estimate at t=k+1');
    xlabel('x'); ylabel('y');

    subplot(3, 2, 4); grid;
    Est_pri2 = mvnpdf([GridX(:) GridY(:)], x_pri(tIdx2)', P_pri(tIdx2,tIdx2));
    Est_pri2 = reshape(Est_pri2, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pri2); grid; title('Apriori Velocity Estimate at t=k+1');
    xlabel('x'); ylabel('y');

    res = y - H*x_pri;
    K = P_pri*H'*(H*P_pri*H'+R)^-1;
    x_pos = x_pri + K*res;
    P_pos = (I_N - K*H)*P_pri;
    subplot(3, 2, 5); 
    Est_pos = mvnpdf([GridX(:) GridY(:)], x_pos(tIdx1)', P_pos(tIdx1,tIdx1));
    Est_pos = reshape(Est_pos, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pos); grid; title('Aposteriori Position Estimate at t=k+1');
    xlabel('x'); ylabel('y');

    subplot(3, 2, 6);
    Est_pos2 = mvnpdf([GridX(:) GridY(:)], x_pos(tIdx2)', P_pos(tIdx2,tIdx2));
    Est_pos2 = reshape(Est_pos2, length(AxisX), length(AxisY));
    contour(AxisX, AxisY, Est_pos2); grid; title('Aposteriori Velocity Estimate at t=k+1');
    xlabel('x'); ylabel('y');
end