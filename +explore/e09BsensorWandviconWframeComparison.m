xsensFName = 'neura-sparse01/calib/S10-Calib-SensorW2V';
viconFName = 'neura-sparse01/calib/S10-Calib-SensorW2V.mat';

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

v_s__q = mocapdb.XsensBody();

w_s__q = mocapdb.XsensBody.loadMTExport(xsensFName, options);
load(viconFName);

m = struct('Pelvis', 'RPV', 'L_UpLeg', 'LTH', 'R_UpLeg', 'RTH', ...
        'L_LowLeg', 'LAK', 'R_LowLeg', 'RAK', ...
        'L_Foot', 'LFT', 'R_Foot', 'RFT');
idx = 20;

for i=1:length(w_s__q.segList)
    n = w_s__q.segList(i); n = n{1};

    if isfield(m, n)
        pN = eval(sprintf('%s_SensN', m.(n)));
        pS = eval(sprintf('%s_SensS', m.(n)));
        pW = eval(sprintf('%s_SensW', m.(n)));

        x = pN(idx, :) -  pS(idx, :);
        x = x / norm(x);

        y = pW(idx, :) - pN(idx, :);
        y = y - dot(x, y)*x;
        y = y / norm(y);

        z = cross(x, y);
        z = z / norm(z);

        v_s__q.(n).ori = rotm2quat([x' y' z']);
    end
end

target = 'Pelvis';
calib = quatmultiply(v_s__q.(target).ori, quatconj(w_s__q.(target).ori(idx, :)));
rad2deg(quat2eul(calib))

target = 'L_LowLeg';
calib = quatmultiply(v_s__q.(target).ori, quatconj(w_s__q.(target).ori(idx, :)));
rad2deg(quat2eul(calib))

updateFigureContents('Rotation'); xlabel('x'); ylabel('y'); zlabel('z'); hold;
pelib.viz.plotR(quat2rotm(w_s__q.(target).ori(idx,:)), [0, 0, 0]);
pelib.viz.plotR(quat2rotm(v_s__q.(target).ori), [0, 0, 0], 'myc');

target = 'R_LowLeg';
calib = quatmultiply(v_s__q.(target).ori, quatconj(w_s__q.(target).ori(idx, :)));
rad2deg(quat2eul(calib))

pelib.viz.plotR(quat2rotm(w_s__q.(target).ori(idx,:)), [-2, 0, 0]);
pelib.viz.plotR(quat2rotm(v_s__q.(target).ori), [-2, 0, 0], 'myc');