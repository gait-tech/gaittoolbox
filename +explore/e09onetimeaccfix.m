%% Fix the faulty hand IMU issue
% S01-Trial
% 90 degrees rotate for left ankle
% {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
%     'JumpingJacks-1', 'JumpingJacks-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
%     'HighKneeJog-1', 'HighKneeJog-2'};
% 270 degrees rotate for left ankle
% {'FigureofEight-1', 'FigureofEight-2', 'Zigzag-1', 'Zigzag-2', ...
%                  'Fivemin-1', 'Fivemin-2'};
% S02-Trial no change

name = 'neura/imu/S02-Trial-';
name_vicon = 'neura/vicon/S02-Trial-';
basis = 'StaticV2W-1';
target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
    'JumpingJacks-1', 'JumpingJacks-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
    'HighKneeJog-1', 'HighKneeJog-2', 'FigureofEight-1', 'FigureofEight-2', ...
    'Zigzag-1', 'Zigzag-2', 'Fivemin-1', 'Fivemin-2'};
% target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
%     'JumpingJacks-1', 'JumpingJacks-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
%     'HighKneeJog-1', 'HighKneeJog-2'};
% target = {'FigureofEight-1', 'FigureofEight-2', 'Zigzag-1', 'Zigzag-2', ...
%                  'Fivemin-1', 'Fivemin-2'};
imus_src = {'00B40C49', '00B40C4A'};
imus_dst = {'00B40C44', '00B40C47'};   
load(sprintf('%s%s', name_vicon, basis));

idx = 5;
quat_viconbasis = {calcQuat(LAK_SensN, LAK_SensS, LAK_SensW, idx), ...
                   calcQuat(RAK_SensN, RAK_SensS, RAK_SensW, idx)};

theta = 0;
R = [cosd(theta) sind(theta) 0;
     -sind(theta) cosd(theta) 0;
     0 0 1];
quat_manual= {rotm2quat(R), rotm2quat(eye(3))};

for i=1:length(imus_src)
    f_basis = sprintf("%s%s-000_%s.txt", name, basis, imus_src{i});
    T_basis = readtable(f_basis);
    quat_basis = [T_basis.Quat_q0(idx) T_basis.Quat_q1(idx) ...
                  T_basis.Quat_q2(idx) T_basis.Quat_q3(idx)];
    
%     quat_adj = quatmultiply(quat_viconbasis{i}, quatconj(quat_basis));
    quat_adj = quat_manual{i};
    for j=1:length(target)
        f_src = sprintf("%s%s-000_%s.txt", name, target{j}, imus_src{i});
        T = readtable(f_src);
        quat_new = [T.Quat_q0(:) T.Quat_q1(:) T.Quat_q2(:) T.Quat_q3(:)];
        quat_new = quatmultiply(quat_adj, quat_new);
        
        T.Quat_q0 = quat_new(:,1);
        T.Quat_q1 = quat_new(:,2);
        T.Quat_q2 = quat_new(:,3);
        T.Quat_q3 = quat_new(:,4);
        
        f_dst = sprintf("%s%s-000_%s.txt", name, target{j}, imus_dst{i});
        writetable(T, f_dst);
    end
end

function q = calcQuat(N, S, W, i)
    x = N(i,:) - S(i,:); x = x / norm(x);
    y = W(i,:) - N(i,:); y = y / norm(y);
    y = y - dot(x, y)*x; y = y / norm(y);
    z = cross(x, y);
    
    q = rotm2quat([x' y' z']);
end