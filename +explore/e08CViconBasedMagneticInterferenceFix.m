global DEGRANGE;
DEGRANGE = (0:0.1:359) - 180;

dir = 'neura-sparse01';
for subjIdx = 1:10
    subj = sprintf('S%02d', subjIdx); 

    options = struct('Pelvis', '00B40B91', ...
        'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
        'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
        'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

    calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
                 sprintf('%s/calib/%s-Calib-SensorW2V.mat', dir, subj), ...
                 sprintf('%s/calib/%s-Calib-SensorW2V', dir, subj), ...
                 options, 100); 
    rad2deg(quat2eul(calibW2V.Pelvis.ori))

    dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Walk-1', dir, subj), options);
    dataS.fs = 100;
    dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Walk-1.mat', dir, subj));
    fs = dataV.fs;

    nSamples = min(dataV.nSamples, dataS.nSamples);
    V__dataV = dataV.getSubset(1:nSamples);
    V__dataV.changePosUnit('m', true);
    W__dataS = dataS.getSubset(1:nSamples);
    % V__dataS = W__dataS.toViconFrame(calibW2V);
    V__dataS = W__dataS;

    sIdx = max(V__dataV.getStartIndex()+1, 100);
    eIdx = length(V__dataV.PELV(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);

    V__viconBody = V__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', V__dataV.fs, ...
                             'xyzColor', {'m', 'y', 'c'}});  
    refAcc = V__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
    refAcc.MIDPEL = refAcc.MIDPEL(sIdx:eIdx,:);
    refAcc.LTIO = refAcc.LTIO(sIdx:eIdx,:);
    refAcc.RTIO = refAcc.RTIO(sIdx:eIdx,:);

    sensors = dataS.exportRawMeasurementAsStruct({'Pelvis', 'L_LowLeg', 'R_LowLeg'}, ...
                        {'PELV', 'LANK', 'RANK'});

    estAcc = {};
    estAcc.MIDPEL = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                            W__dataS.Pelvis.acc) - [0 0 9.81];
    estAcc.MIDPEL = estAcc.MIDPEL(sIdx:eIdx,:);
    % estAcc.MIDPEL = quatrotate(quatconj(calibW2V.Pelvis.ori), estAcc.MIDPEL);
    estAcc.LTIO = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                            W__dataS.L_LowLeg.acc) - [0 0 9.81];
    estAcc.LTIO = estAcc.LTIO(sIdx:eIdx,:);
    % estAcc.LTIO = quatrotate(quatconj(calibW2V.L_LowLeg.ori), estAcc.LTIO);
    estAcc.RTIO = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                            W__dataS.R_LowLeg.acc) - [0 0 9.81];
    estAcc.RTIO = estAcc.RTIO(sIdx:eIdx,:);
    % estAcc.RTIO = quatrotate(quatconj(calibW2V.R_LowLeg.ori), estAcc.RTIO);

    calibThetaL = rad2deg(quat2eul(quatconj(calibW2V.L_LowLeg.ori)));
    calibThetaR = rad2deg(quat2eul(quatconj(calibW2V.R_LowLeg.ori)));

    calib2 = dataS.calcCalibAnkleSensorW2PelvisWFromAcc(1:1300);
    calib2ThetaL = rad2deg(quat2eul(quatconj(calib2.L_LowLeg.ori)));
    calib2ThetaR = rad2deg(quat2eul(quatconj(calib2.R_LowLeg.ori)));

    viconCalibSB = V__dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
    calib3 = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
    calib3ThetaL = rad2deg(quat2eul(quatconj(calib3.L_LowLeg.ori)));
    calib3ThetaR = rad2deg(quat2eul(quatconj(calib3.R_LowLeg.ori)));

    calib4 = dataS.calcCalibAnkleSensorW2PelvisWFromGyroSkewness(DEGRANGE);
    calib4ThetaL = rad2deg(quat2eul(quatconj(calib4.L_LowLeg.ori)));
    calib4ThetaR = rad2deg(quat2eul(quatconj(calib4.R_LowLeg.ori)));

    [thetaL, errL] = findOptimalThetaBrute(refAcc.LTIO, estAcc.LTIO);
    [thetaR, errR] = findOptimalThetaBrute(refAcc.RTIO, estAcc.RTIO);

    updateFigureContents('Finding Theta'); hold on;
    [~, peakIdxL] = min(abs(DEGRANGE-thetaL)); 
    [~, calibIdxL] = min(abs(DEGRANGE-calibThetaL(1)));
    [~, calib2IdxL] = min(abs(DEGRANGE-calib2ThetaL(1)));
    [~, calib3IdxL] = min(abs(DEGRANGE-calib3ThetaL(1)));
    [~, calib4IdxL] = min(abs(DEGRANGE-calib4ThetaL(1)));
    plot(DEGRANGE, errL, 'r-x', 'MarkerIndices', [peakIdxL]);
    plot(calibThetaL(1), errL(calibIdxL), 'r+');
    plot(calib2ThetaL(1), errL(calib2IdxL), 'ro');
    plot(calib3ThetaL(1), errL(calib3IdxL), 'r*');
    plot(calib4ThetaL(1), errL(calib4IdxL), 'r^');
    [~, peakIdxR] = min(abs(DEGRANGE-thetaR)); 
    [~, calibIdxR] = min(abs(DEGRANGE-calibThetaR(1)));
    [~, calib2IdxR] = min(abs(DEGRANGE-calib2ThetaR(1)));
    [~, calib3IdxR] = min(abs(DEGRANGE-calib3ThetaR(1)));
    [~, calib4IdxR] = min(abs(DEGRANGE-calib4ThetaR(1)));
    plot(DEGRANGE, errR, 'b-x', 'MarkerIndices', [peakIdxR]);
    plot(calibThetaR(1), errR(calibIdxR), 'b+');
    plot(calib2ThetaR(1), errR(calib2IdxR), 'bo');
    plot(calib3ThetaR(1), errR(calib3IdxR), 'b*');
    plot(calib4ThetaR(1), errR(calib4IdxR), 'b^');

    % legend(sprintf('Left  Opt at %.1f%c (%.2f)', thetaL, char(176), errL(peakIdxL)), ...
    %        sprintf('Right Opt at %.1f%c (%.2f)', thetaR, char(176), errR(peakIdxR)) );
    legend(sprintf('Left  vicon at %.1f%c (%.2f)', thetaL, char(176), errL(peakIdxL)), ...
           sprintf('Left  tracking mount at %.1f%c (%.2f)', calibThetaL(1), char(176), errL(calibIdxL)), ...
           sprintf('Left  walk (acc) at %.1f%c (%.2f)', calib2ThetaL(1), char(176), errL(calib2IdxL)), ...
           sprintf('Left  walk (ROM) at %.1f%c (%.2f)', calib3ThetaL(1), char(176), errL(calib3IdxL)), ...
           sprintf('Left  walk (gyro skewness) at %.1f%c (%.2f)', calib4ThetaL(1), char(176), errL(calib4IdxL)), ...
           sprintf('Right vicon at %.1f%c (%.2f)', thetaR, char(176), errR(peakIdxR)), ...
           sprintf('Right tracking mount at %.1f%c (%.2f)', calibThetaR(1), char(176), errR(calibIdxR)), ...
           sprintf('Right walk (acc) at %.1f%c (%.2f)', calib2ThetaR(1), char(176), errR(calib2IdxR)), ...
           sprintf('Right walk (ROM) at %.1f%c (%.2f)', calib3ThetaR(1), char(176), errR(calib3IdxR)), ...
           sprintf('Right walk (gyro skewness) at %.1f%c (%.2f)', calib4ThetaR(1), char(176), errR(calib4IdxR)) );
    title(sprintf('%s', subj));
    
    saveas(gcf, sprintf('explore_output/MagneticInterferenceFix-%s.png', subj))
end

% updateFigureContents('Left  Free Body Acc (World Frame)');
% est = estAcc.LTIO; ref = refAcc.LTIO; 
% estOpt = rotateThetaAlongZ(estAcc.LTIO, thetaL);
% estCal = rotateThetaAlongZ(estAcc.LTIO, calibThetaL(1));
% pelib.viz.plotXYZ(fs, ref, est, estOpt, estCal);
% 
% updateFigureContents('Right Free Body Acc (World Frame)');
% est = estAcc.RTIO; ref = refAcc.RTIO; 
% estOpt = rotateThetaAlongZ(estAcc.RTIO, thetaR);
% estCal = rotateThetaAlongZ(estAcc.RTIO, calibThetaR(1));
% pelib.viz.plotXYZ(fs, ref, est, estOpt, estCal);

% lqdiff = zeros(size(W__dataS.L_LowLeg.ori, 1)-1, 4);
% for i=1:size(lqdiff, 1)
%     lqdiff(i, :) = quatmultiply(W__dataS.L_LowLeg.ori(i+1, :), ...
%                         quatconj(W__dataS.L_LowLeg.ori(i, :)));
% end
% lqdiff = lqdiff((sIdx:eIdx)-1, :);
% lqdiff2 = quat2axang(lqdiff);

% lqdiff = sensors.LANKGyr(sIdx:eIdx, :)
% VAR_WIN  = floor(fs*0.25);
% movVarAcc = movingvar(sqrt( sum(estAcc.LTIO .^2, 2)), VAR_WIN);
% bIsStep = movVarAcc < 1;
% lqdiff2 = lqdiff;
% lqdiff2(bIsStep, :) = NaN;
% plot(lqdiff2)

function [theta, err] = findOptimalThetaBrute(refAcc, estAcc)
    global DEGRANGE;
    err = zeros(length(DEGRANGE), 1);
    for i=1:length(DEGRANGE)
        estAcc2 = rotateThetaAlongZ(estAcc, DEGRANGE(i));
        err(i) = sqrt(mean(nanmean((refAcc(:,1:2) - estAcc2(:,1:2)).^2, 1)));
    end
    [M, I] = min(err);
    theta = DEGRANGE(I);
end

function newData = rotateThetaAlongZ(data, deg)
    R = [cosd(deg) sind(deg) 0;
         -sind(deg) cosd(deg) 0;
             0 0 1];
    newData = quatrotate(quatconj(rotm2quat(R)), data);
end