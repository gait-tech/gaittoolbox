
dir = 'neura-sparse01';
subj = 'S10';

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
             sprintf('%s/calib/%s-Calib-SensorW2V.mat', dir, subj), ...
             sprintf('%s/calib/%s-Calib-SensorW2V', dir, subj), ...
             options, 100); 
dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/calib/%s-Calib-SensorW2V', dir, subj), options);

LAKMag = quatrotate(quatconj(dataS.L_LowLeg.ori), dataS.L_LowLeg.mag);
RAKMag = quatrotate(quatconj(dataS.R_LowLeg.ori), dataS.R_LowLeg.mag);
RPVMag = quatrotate(quatconj(dataS.Pelvis.ori), dataS.Pelvis.mag);

theta1 = acosd(dot(LAKMag(1,1:2), RPVMag(1,1:2))/(norm(LAKMag(1,1:2))*norm(RPVMag(1,1:2))));
R1 = [cosd(theta1) sind(theta1) 0;
     -sind(theta1) cosd(theta1) 0;
     0 0 1];
theta2 = acosd(dot(RAKMag(1,1:2), RPVMag(1,1:2))/(norm(RAKMag(1,1:2))*norm(RPVMag(1,1:2))));
R2 = [cosd(theta2) sind(theta2) 0;
     -sind(theta2) cosd(theta2) 0;
     0 0 1];
[theta1 theta2]

LAKOri2 = quatmultiply(rotm2quat(R1), dataS.L_LowLeg.ori);
RAKOri2 = quatmultiply(rotm2quat(R2), dataS.R_LowLeg.ori);
% LAKOri2 = quatmultiply(calibW2V.L_LowLeg.ori, dataS.L_LowLeg.ori);
% RAKOri2 = quatmultiply(calibW2V.R_LowLeg.ori, dataS.R_LowLeg.ori);
LAKMag2 = quatrotate(quatconj(LAKOri2), dataS.L_LowLeg.mag);
RAKMag2 = quatrotate(quatconj(RAKOri2), dataS.R_LowLeg.mag);

updateFigureContents('Magnetometer (World Frame)');
pelib.viz.plotXYZ(100, RPVMag, LAKMag, RAKMag, LAKMag2, RAKMag2);


