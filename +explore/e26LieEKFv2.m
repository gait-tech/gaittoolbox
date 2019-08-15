R0 = eul2rotm([pi/2 0 pi/2]);
p0 = [0 1 1]';
T0 = [R0 p0; zeros(1,3) 1];

R1 = eul2rotm([pi/100 pi/100 0]);
p1 = [0 0.5 1]';
T1 = [R1 p1; zeros(1,3) 1];

T2a = T0*T1
T2b = T1*T0