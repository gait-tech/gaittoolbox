name = 'xsensAnkleExperiment\MT_012000EB-005';
options = struct('R_UpLeg', '00B40C3C', 'R_LowLeg', '00B40C35', 'R_Ankle', '00B40C48');
bs = fieldnames(options);
bsN = length(bs);

data = {};
for i = 1:bsN
    fpath = sprintf("%s-000_%s.txt", name, options.(bs{i}));
    T = readtable(fpath);
    data.(bs{i}) = table([T.Quat_q0 T.Quat_q1 T.Quat_q2 T.Quat_q3], ...
                [T.Acc_X T.Acc_Y T.Acc_Z], ...
                [T.Gyr_X T.Gyr_Y T.Gyr_Z], ...
                [T.Mag_X T.Mag_Y T.Mag_Z], ...
               'VariableNames', {'ori', 'acc', 'gyr', 'mag'});
end
qI = rotm2quat(eye(3));

out = {};
for i = 1:bsN
   S_q_B = quatmultiply(quatconj(data.(bs{i}).ori(1,:)), qI);
   out.(bs{i}) = quatmultiply(data.(bs{i}).ori, S_q_B);
end

RUpLeg = out.R_UpLeg; RLowLeg = out.R_LowLeg; RAnkle = out.R_Ankle;
updateFigureContents('RLeg Quaternion');
pelib.viz.plotQuaternion(RUpLeg, RLowLeg, RAnkle);

RKneeAngle_RLowLeg = pelib.grBody.calcJointAngles(RUpLeg, RLowLeg)*180/pi;
RKneeAngle_RAnkle = pelib.grBody.calcJointAngles(RUpLeg, RAnkle)*180/pi;
updateFigureContents('RKnee Joint Angles');
pelib.viz.plotXYZ(RKneeAngle_RLowLeg, RKneeAngle_RAnkle)