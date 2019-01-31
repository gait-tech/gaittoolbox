% motion list
listSetup = struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M51+C373");

dataSfname = sprintf('neura-sparse01/imu/%s', listSetup.file);
load(sprintf('neura-sparse01/explore-v2/neura-%s-debug.mat', listSetup.file));
load(sprintf('neura-sparse01/explore-v2/neura-%s-%s.mat', listSetup.file, listSetup.algo));
    
targetname = sprintf('explore_output/neura-%s-%s', listSetup.file, listSetup.algo);

if cs.initSrc == 'w__v'
    aLabel = 'w__v';
    vb = W__viconBody;
elseif cs.initSrc == 'v__v'
    aLabel = 'v__v';
    vb = V__viconBody;
else
    aLabel = 'w__x';
    vb = W__xsensBody;
end
if ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
   strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
    eLabel = strcat(cs.accData, cs.initSrc(end));
else
    eLabel = cs.accData;
end
idx = allIdx.(cs.initSrc);

%% Initialize estState, estP, estBody, and actBody debug versions
n2 = size(estState, 1)*3;
estStateDebug = repelem(estState, 3, 1);
estStateDebug(1:3:end, :) = estState2.predState;
estStateDebug(2:3:end, :) = estState2.zuptState;
estPDebug = repelem(estState2.predP, 1, 1, 3);
estPDebug(:, :, 2:3:end) = estState2.zuptP;
estPDebug(:, :, 3:3:end) = estState2.cstrP;

% options = struct('Pelvis', '00B40B91', ...
%     'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%     'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
%     'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
% dataS = mocapdb.XsensBody.loadMTExport(dataSfname, options);
% dataS = dataS.getSubset(idx);

estBodyDebug = estBody.copy();
estBodyDebug.MIDPEL = repelem(estBodyDebug.MIDPEL, 3, 1);
estBodyDebug.MIDPEL(1:3:end, :) = estState2.predState(:,1:3);
estBodyDebug.MIDPEL(2:3:end, :) = estState2.zuptState(:,1:3);
estBodyDebug.LTIO = repelem(estBodyDebug.LTIO, 3, 1);
estBodyDebug.LTIO(1:3:end, :) = estState2.predState(:,11:13);
estBodyDebug.LTIO(2:3:end, :) = estState2.zuptState(:,11:13);
estBodyDebug.RTIO = repelem(estBodyDebug.RTIO, 3, 1);
estBodyDebug.RTIO(1:3:end, :) = estState2.predState(:,21:23);
estBodyDebug.RTIO(2:3:end, :) = estState2.zuptState(:,21:23);

estBodyDebug.qRPV = repelem(estBodyDebug.qRPV, 3, 1);
estBodyDebug.qLSK = repelem(estBodyDebug.qLSK, 3, 1);
estBodyDebug.qRSK = repelem(estBodyDebug.qRSK, 3, 1);

v = quat2rotm(estBodyDebug.qLSK); v = squeeze(v(:,3,:))';
estBodyDebug.LFEO = estBodyDebug.LTIO + estBodyDebug.getLShankLength()*v;
v = quat2rotm(estBodyDebug.qRSK); v = squeeze(v(:,3,:))';
estBodyDebug.RFEO = estBodyDebug.RTIO + estBodyDebug.getRShankLength()*v;
v = quat2rotm(estBodyDebug.qRPV); v = squeeze(v(:,2,:))';
estBodyDebug.LFEP = estBodyDebug.MIDPEL + estBodyDebug.getPelvisLength()/2*v;
estBodyDebug.RFEP = estBodyDebug.MIDPEL - estBodyDebug.getPelvisLength()/2*v;

v = zeros(3, 3, n2);
z = (estBodyDebug.LFEP-estBodyDebug.LFEO)';
z = z ./ vecnorm(z, 2, 1);
v(:, 3, :) = reshape(z, 3, 1, []);
y =  quat2rotm(estBodyDebug.qLSK);
v(:, 2, :) = y(:, 2, :);
x = cross(v(:, 2, :), v(:, 3, :));
x = x ./ vecnorm(x, 2, 1);
v(:, 1, :) =  reshape(x, 3, 1, []);
estBodyDebug.qLTH = rotm2quat(v);

v = zeros(3, 3, n2);
z = (estBodyDebug.RFEP-estBodyDebug.RFEO)';
z = z ./ vecnorm(z, 2, 1);
v(:, 3, :) = reshape(z, 3, 1, []);
y =  quat2rotm(estBodyDebug.qRSK);
v(:, 2, :) = y(:, 2, :);
x = cross(v(:, 2, :), v(:, 3, :));
x = x ./ vecnorm(x, 2, 1);
v(:, 1, :) =  reshape(x, 3, 1, []);
estBodyDebug.qRTH = rotm2quat(v);

estBodyDebug.nSamples = estBodyDebug.nSamples*3;
estBodyDebug.fs = estBodyDebug.fs*3;

viconBodyDebug = vb.copy();
for i=1:length(viconBodyDebug.posList)
    n = viconBodyDebug.posList{i};
    viconBodyDebug.(n) = repelem(viconBodyDebug.(n), 3, 1);
end
for i=1:length(viconBodyDebug.oriList)
    n = viconBodyDebug.oriList{i};
    viconBodyDebug.(n) = repelem(viconBodyDebug.(n), 3, 1);
end
viconBodyDebug.nSamples = viconBodyDebug.nSamples*3;
viconBodyDebug.LTOE = [];
viconBodyDebug.RTOE = [];

estPMPPosDebug = estPDebug( 1: 3,  1: 3, :);
estPLAPosDebug = estPDebug(11:13, 11:13, :);
estPRAPosDebug = estPDebug(21:23, 21:23, :);

%% Summary Values
% plot P and rank
idx2 = 1:estBodyDebug.nSamples;
estPRank = zeros(estBodyDebug.nSamples, 1);
for i=1:n2
    estPRank(i,1) = rank(estPDebug(:,:,i));
end

updateFigureContents('P and rank');
subplot(4, 1, 1);
plot(idx2, squeeze(estPMPPosDebug(1,1,:)), 'r', ...
     idx2, squeeze(estPMPPosDebug(2,2,:)), 'g', ...
     idx2, squeeze(estPMPPosDebug(3,3,:)), 'b');
title('MP P diagonals');
legend('x', 'y', 'z');

subplot(4, 1, 2);
plot(idx2, squeeze(estPLAPosDebug(1,1,:)), 'r', ...
     idx2, squeeze(estPLAPosDebug(2,2,:)), 'g', ...
     idx2, squeeze(estPLAPosDebug(3,3,:)), 'b');
title('LA P diagonals');
legend('x', 'y', 'z');

subplot(4, 1, 3);
plot(idx2, squeeze(estPRAPosDebug(1,1,:)), 'r', ...
     idx2, squeeze(estPRAPosDebug(2,2,:)), 'g', ...
     idx2, squeeze(estPRAPosDebug(3,3,:)), 'b');
title('RA P diagonals');
legend('x', 'y', 'z');

subplot(4, 1, 4);
plot(idx2, estPRank);
title('Rank');

% plot rank and femur length
updateFigureContents('Femur length'); 
subplot(2, 1, 1);
LFemurLength = estBodyDebug.calcLFemurLength();
plot(idx2, LFemurLength, ...
     idx2, repelem(LFemurLength(1), length(idx2)));
legend('est LFem len', 'act LFem len');
RFemurLength = estBodyDebug.calcRFemurLength();
subplot(2, 1, 2);
plot(idx2, RFemurLength, ...
     idx2, repelem(RFemurLength(1), length(idx2)));
legend('est RFem len', 'act RFem len');

vel1 = vb.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
vel2 = estBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});

updateFigureContents('Vicon Joint Velocities'); 
viconMPVel = vel1.MIDPEL; viconLAVel = vel1.LTIO; viconRAVel = vel1.RTIO;
viconLARAMeanVel = (viconLAVel+viconRAVel)./2;
pelib.viz.plotXYZ(fs, viconMPVel, viconLARAMeanVel, viconLAVel, viconRAVel);

updateFigureContents('Vicon Joint Positions'); 
viconMPPos = vb.MIDPEL; viconLAPos = vb.LTIO; viconRAPos = vb.RTIO;
viconLARAMeanPos = (viconLAPos+viconRAPos)./2;
pelib.viz.plotXYZ(fs, viconMPPos, viconLARAMeanPos, viconLAPos, viconRAPos);

%% Animation
% az = 0; el = 180;
updateFigureContents('Animation'); 
tmpBody1 = estBodyDebug;
tmpBody2 = false;
tmpBody1Limits = [tmpBody1.xlim() tmpBody1.ylim() tmpBody1.zlim()];
samplePlotN = 1000;
for i=1:2844
    [az, el] = view;
    clf; grid; axis square;
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim(tmpBody1Limits(1:2)); 
    ylim(tmpBody1Limits(3:4)); 
    zlim(tmpBody1Limits(5:6));  
    view(az, el);
    pelib.viz.plotLowerBody(tmpBody1, i, true, false);
    hold on;
    pelib.viz.plotMeanCovSamples(tmpBody1.RTIO(i, :), ...
                estPRAPosDebug(:, :, i), samplePlotN, ...
                sprintf('%s.', tmpBody1.xyzColor{1}));
    pelib.viz.plotMeanCovSamples(tmpBody1.MIDPEL(i, :), ...
                estPMPPosDebug(:, :, i), samplePlotN, ...
                sprintf('%s.', tmpBody1.xyzColor{2}));
    pelib.viz.plotMeanCovSamples(tmpBody1.LTIO(i, :), ...
                estPLAPosDebug(:, :, i), samplePlotN, ...
                sprintf('%s.', tmpBody1.xyzColor{3}));
    if ~(tmpBody2==false)
        pelib.viz.plotLowerBody(tmpBody2, i, true, false);
    end
    
    pause(1/1000);
end