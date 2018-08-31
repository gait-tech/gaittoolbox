%% Acc Exploration 2 (neura dataset)
fs=100;
optionTIB = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

% optionANK = struct('Pelvis', '00B40B91', ...
%     'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%     'L_LowLeg', '00B40C49', 'R_LowLeg', '00B40C4A', ...
%     'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
optionANK = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
fnameV = 'neura/vicon/S03-Trial-000.mat';
fnameX = 'neura/xsens/S01-Trial-006.bvh';
fnameS = 'neura/imu/S03-Trial-000';

% Rw2v = mocapdb.loadPendulumCompassMat('neura\calib\Pendulum02.mat', 'neura\calib\XsensCompass02.mat');
viconANK = mocapdb.ViconBody.loadViconMat(fnameV);
xsensANK = mocapdb.BVHBody.loadXsensBVHFile(fnameX, "mm");

% viconANK = viconANK.toWorldFrame(rotm2quat(Rw2v'));
sensTIB = mocapdb.XsensBody.loadMTExport(fnameS, optionTIB);
sensANK = mocapdb.XsensBody.loadMTExport(fnameS, optionANK);

idx = 1:(viconANK.nSamples-1);
% viconBody = viconANK.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
%                              'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
%                              'xyzColor', {'m', 'y', 'c'}}); 
% viconBodyRel = viconBody.changeRefFrame('MIDPEL');
% viconCalibWB = sensANK.calcCalibSB(viconBody, sIdx);
    
gfrAcc = {};
% viconANK
val1 = {'PELV', 'LFEP', 'RFEP', 'LFEO', 'RFEO', 'LTIO', 'RTIO', 'RTOE', 'LTOE'};
for i=1:length(val1)
    viconANK.(val1{i}) = viconANK.(val1{i})./1000;
end
MIDPEL_vicon = [mean([viconANK.LFEP(:,1) viconANK.RFEP(:,1)], 2),...
                mean([viconANK.LFEP(:,2) viconANK.RFEP(:,2)], 2),...
                mean([viconANK.LFEP(:,3) viconANK.RFEP(:,3)], 2)];
sIdx = 5; eIdx = min(length(MIDPEL_vicon(:,1))-1, sensTIB.nSamples-1);
gfrAccvANKMP = [0 0 0; diff(MIDPEL_vicon, 2, 1)*fs*fs];
gfrAccvANKMP = gfrAccvANKMP(sIdx:eIdx,:);
gfrAccvANKLA = [0 0 0; diff(viconANK.LTIO, 2, 1)*fs*fs];
gfrAccvANKLA = gfrAccvANKLA(sIdx:eIdx,:);
gfrAccvANKRA = [0 0 0; diff(viconANK.RTIO, 2, 1)*fs*fs];
gfrAccvANKRA = gfrAccvANKRA(sIdx:eIdx,:);              

% xsensANK
val1 = {'Hips', 'RightUpLeg', 'RightLeg', 'RightFoot', 'RightToe', ...
        'LeftUpLeg', 'LeftLeg', 'LeftFoot', 'LeftToe'};
for i=1:length(val1)
    xsensANK.(val1{i}) = xsensANK.(val1{i})./1000;
end
MIDPEL_xsens = [mean([xsensANK.LeftUpLeg(:,1) xsensANK.RightUpLeg(:,1)], 2),...
                mean([xsensANK.LeftUpLeg(:,2) xsensANK.RightUpLeg(:,2)], 2),...
                mean([xsensANK.LeftUpLeg(:,3) xsensANK.RightUpLeg(:,3)], 2)];
gfrAccxANKMP = [0 0 0; diff(MIDPEL_xsens, 2, 1)*fs*fs];
gfrAccxANKMP = gfrAccxANKMP(sIdx:eIdx,:);
gfrAccxANKLA = [0 0 0; diff(xsensANK.LeftFoot, 2, 1)*fs*fs];
gfrAccxANKLA = gfrAccxANKLA(sIdx:eIdx,:);
gfrAccxANKRA = [0 0 0; diff(xsensANK.RightFoot, 2, 1)*fs*fs];
gfrAccxANKRA = gfrAccxANKRA(sIdx:eIdx,:); 
        
% sensTIB
gfrAccsTIBMP = quatrotate(quatconj(sensTIB.Pelvis.ori), ...
                            sensTIB.Pelvis.acc) - [0 0 9.81];
gfrAccsTIBMP = gfrAccsTIBMP(sIdx:eIdx,:);
gfrAccsTIBLA = quatrotate(quatconj(sensTIB.L_LowLeg.ori), ...
                            sensTIB.L_LowLeg.acc) - [0 0 9.81];
gfrAccsTIBLA = gfrAccsTIBLA(sIdx:eIdx,:);
gfrAccsTIBRA = quatrotate(quatconj(sensTIB.R_LowLeg.ori), ...
                            sensTIB.R_LowLeg.acc) - [0 0 9.81];
gfrAccsTIBRA = gfrAccsTIBRA(sIdx:eIdx,:);

% sensANK
% qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
% xsensANK.Pelvis.ori = quatmultiply(xsensANK.Pelvis.ori, qTCD2BM);
% xsensANK.L_LowLeg.ori = quatmultiply(xsensANK.L_LowLeg.ori, qTCD2BM);
% xsensANK.R_LowLeg.ori = quatmultiply(xsensANK.R_LowLeg.ori, qTCD2BM);

% theta = 180;
% R = [cosd(theta) sind(theta) 0;
%      -sind(theta) cosd(theta) 0;
%      0 0 1];
% sensANK.L_LowLeg.ori = quatmultiply(rotm2quat(R), sensANK.L_LowLeg.ori);
% theta = 90;
% R = [cosd(theta) sind(theta) 0;
%      -sind(theta) cosd(theta) 0;
%      0 0 1];
% sensANK.R_LowLeg.ori = quatmultiply(rotm2quat(R), sensANK.R_LowLeg.ori);
  
gfrAccsANKMP = quatrotate(quatconj(sensANK.Pelvis.ori), ...
                            sensANK.Pelvis.acc) - [0 0 9.81];
gfrAccsANKMP = gfrAccsANKMP(sIdx:eIdx,:);
gfrAccsANKLA = quatrotate(quatconj(sensANK.L_LowLeg.ori), ...
                            sensANK.L_LowLeg.acc) - [0 0 9.81];
gfrAccsANKLA = gfrAccsANKLA(sIdx:eIdx,:);
gfrAccsANKRA = quatrotate(quatconj(sensANK.R_LowLeg.ori), ...
                            sensANK.R_LowLeg.acc) - [0 0 9.81];
gfrAccsANKRA = gfrAccsANKRA(sIdx:eIdx,:);

updateFigureContents('Hip Acc');
% clf; pelib.viz.plotXYZ(fs, gfrAccvANKMP, gfrAccxANKMP, gfrAccsANKMP, gfrAccsTIBMP);
clf; pelib.viz.plotXYZ(fs, gfrAccvANKMP, gfrAccsANKMP, gfrAccsTIBMP);

updateFigureContents('LA Acc');
% clf; pelib.viz.plotXYZ(fs, gfrAccvANKLA, gfrAccxANKLA, gfrAccsANKLA, gfrAccsTIBLA);
clf; pelib.viz.plotXYZ(fs, gfrAccvANKLA, gfrAccsANKLA, gfrAccsTIBLA);

updateFigureContents('RA Acc');
% clf; pelib.viz.plotXYZ(fs, gfrAccvANKRA, gfrAccxANKRA, gfrAccsANKRA, gfrAccsTIBRA);
clf; pelib.viz.plotXYZ(fs, gfrAccvANKRA, gfrAccsANKRA, gfrAccsTIBRA);

az = 0; el = 180;
updateFigureContents('Animation');
sIdx = 20;
tmpBody1 = viconANK.togrBody(sIdx:sIdx+1, {});
tmpBody1.qRPV = sensANK.Pelvis.ori(sIdx:end,:);
tmpBody1.qLTH = sensANK.L_UpLeg.ori(sIdx:end,:);
tmpBody1.qRTH = sensANK.R_UpLeg.ori(sIdx:end,:);
tmpBody1.qLSK = sensANK.L_LowLeg.ori(sIdx:end,:);
tmpBody1.qRSK = sensANK.R_LowLeg.ori(sIdx:end,:);
tmpBody1.qLFT = sensANK.L_Foot.ori(sIdx:end,:);
tmpBody1.qRFT = sensANK.R_Foot.ori(sIdx:end,:);
tmpBody2 = false;
tmpBody1Limits = [tmpBody1.xlim()+[-1 1] tmpBody1.ylim()+[-1 1] tmpBody1.zlim()+[-1 1]];
i = 1; 
while i <= tmpBody1.nSamples
    [az, el] = view;
    clf; grid;
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim(tmpBody1Limits(1:2)); 
    ylim(tmpBody1Limits(3:4)); 
    zlim(tmpBody1Limits(5:6));  
    view(az, el);
    pelib.viz.plotLowerBody(tmpBody1, i, true, false);
    if ~(tmpBody2==false)
        pelib.viz.plotLowerBody(tmpBody2, i, true, false);
    end
    i = i+5;
    pause(1/1000);
end