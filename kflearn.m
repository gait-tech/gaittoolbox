dt = 1;
I_N = eye(4);

F = [1 0 dt 0; 0 0 1 0; 0 1 0 dt; 0 0 0 1];
G = [0.5*dt^2 0; 0 0.5*dt^2; dt 0; 0 dt];
H = [0 1 0 0];

Q = [0.01 0 0 0; 0 0.01 0 0; 0 0 0.01 0; 0 0 0 0.01];
R = [1e-1];

x_pos = [0; 0; 0; 0];
u = [1; 1];
P_pos = [1 0 0 0; 0 1 0 0; 0 0 2 0; 0 0 0 2];

x_pri = F*x_pos+G*u;
P_pri = F * P_pos * F' + Q;

res = [1] - H*x_pri;
K = P_pri*H'*(H*P_pri*H'+R)^-1;
x_pos = x_pri + K*res;
P_pos = (I_N - K*H)*P_pri;