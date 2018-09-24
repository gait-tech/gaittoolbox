% Show the absolute and relative position
fs=100;
n = viconBody.nSamples;
actState = [viconBody.MIDPEL zeros(n, 3) qOri.v.PELV...
           viconBody.LTIO zeros(n, 3) qOri.v.LTIB...
           viconBody.RTIO zeros(n, 3) qOri.v.RTIB];
               
updateFigureContents('Hips');
clf; pelib.viz.plotStateComparison(estState2, actState, 1, fs);

updateFigureContents('LAnkle');
clf; pelib.viz.plotStateComparison(estState2, actState, 11, fs);

updateFigureContents('RAnkle');
clf; pelib.viz.plotStateComparison(estState2, actState, 21, fs);

updateFigureContents('LSHK Quaternion');
clf; pelib.viz.plotOrientation({viconBody, estBody}, {'qLSK', 'qRSK'})

% estStateRel2 = pelib.est.changeStateRefFrame(estState2);
% actStateRel = pelib.est.changeStateRefFrame(actState);
% updateFigureContents('LAnkleRel');
% clf; pelib.viz.plotStateComparison(estStateRel2, actStateRel, 11, 60);
% 
% updateFigureContents('RAnkleRel');
% clf; pelib.viz.plotStateComparison(estStateRel2, actStateRel, 21, 60);
% 
% updateFigureContents('Hips2');
% clf; pelib.viz.plotStateComparison(estState2, actState, 1, 60);
% 
% updateFigureContents('LAnkle2');
% clf; pelib.viz.plotStateComparison(estState2, actState, 11, 60);
% 
% updateFigureContents('RAnkle2');
% clf; pelib.viz.plotStateComparison(estState2, actState, 21, 60);
% 
% estStateRel2 = pelib.est.changeStateRefFrame(estState2);
% actStateRel = pelib.est.changeStateRefFrame(actState);
% updateFigureContents('LAnkleRel2');
% clf; pelib.viz.plotStateComparison(estStateRel2, actStateRel, 11, 60);
% 
% updateFigureContents('RAnkleRel2');
% clf; pelib.viz.plotStateComparison(estStateRel2, actStateRel, 21, 60);
% 
% % P exploration
% load('C:\Users\z5151460\workspace\gaitrecon\experiments\tcd-s1-acting1-Dxxx+ZUPT+C001.mat')
% [nSamples, nState] = size(estState);
% 
% P_pri1 = nan(nSamples, nState);
% P_pos1 = nan(nSamples, nState);
% P_con1 = nan(nSamples, nState);
% 
% for i=1:nSamples
%     P_pri1(i,:) = diag(estState2.predP(:,:,i));
%     P_pos1(i,:) = diag(estState2.zuptP(:,:,i));
%     P_con1(i,:) = diag(estState2.cstrP(:,:,i));
% end
% 
% load('C:\Users\z5151460\workspace\gaitrecon\experiments\tcd-s1-acting1-Dxxx+ZUPT+C002.mat')
% [nSamples, nState] = size(estState);
% 
% P_pri2 = nan(nSamples, nState);
% P_pos2 = nan(nSamples, nState);
% P_con2 = nan(nSamples, nState);
% 
% for i=1:nSamples
%     P_pri2(i,:) = diag(estState2.predP(:,:,i));
%     P_pos2(i,:) = diag(estState2.zuptP(:,:,i));
%     P_con2(i,:) = diag(estState2.cstrP(:,:,i));
% end

%% Static Plots
% preprocessing
viconBodyRel = viconBody.changeRefFrame('MIDPEL');
estBodyRel = estBody.changeRefFrame('MIDPEL');
xsensBodyRel = xsensBody.changeRefFrame('MIDPEL');

updateFigureContents('Position');
pelib.viz.plotPosition({viconBodyRel, xsensBodyRel}, {'LTIO', 'RTIO'});

updateFigureContents('Position');
pelib.viz.plotPosition({viconBodyRel, estBodyRel}, {'LTIO', 'RTIO'});


viconBodyRel = viconBody.changeRefFrame('MIDPEL');
viconBodyRel.xyzColor = {'k', 'k', 'k'};
load('C:\Users\z5151460\workspace\gaitrecon\neura\explore\neura-S03-Trial-Fivemin-1-Nssv+M02+C202.mat')
estBodyRel_Lin = estBody.changeRefFrame('MIDPEL');
estBodyRel_Lin.name = 'lin est'; estBodyRel_Lin.xyzColor = {'r','r','r'};
load('C:\Users\z5151460\workspace\gaitrecon\neura\explore\neura-S03-Trial-Fivemin-1-Nssv+M02+C152.mat')
estBodyRel_Nonlin = estBody.changeRefFrame('MIDPEL'); 
estBodyRel_Nonlin.name = 'nonlin est'; estBodyRel_Nonlin.xyzColor = {'b','b','b'};

updateFigureContents('Position');
pelib.viz.plotPosition({viconBodyRel, estBodyRel_Lin, estBodyRel_Nonlin}, {'LTIO', 'RTIO'});

updateFigureContents('Joint Angles (Knee)');
pelib.viz.plotJointAngles({viconBodyRel, estBodyRel_Lin, estBodyRel_Nonlin}, {'LKnee', 'RKnee'})

updateFigureContents('Joint Angles (Hips)');
pelib.viz.plotJointAngles({viconBodyRel, xsensBodyRel}, {'LHip', 'RHip'})

updateFigureContents('Joint Angles (Hips)');
pelib.viz.plotJointAngles({viconBodyRel, estBodyRel}, {'LHip', 'RHip'})

updateFigureContents('Joint Angles (Knee)');
pelib.viz.plotJointAngles({viconBodyRel, estBodyRel}, {'LKnee', 'RKnee'})

updateFigureContents('Joint Angles (Knee)');
pelib.viz.plotJointAngles({viconBodyRel, xsensBodyRel}, {'LKnee', 'RKnee'})

% Constraint Info
% d_pelvis = norm(actBody.LFEP(1,:) - actBody.RFEP(1,:));
% d_lfemur = norm(actBody.LFEP(1,:) - actBody.LFEO(1,:));
% d_rfemur = norm(actBody.RFEP(1,:) - actBody.RFEO(1,:));
% d_ltibia = norm(actBody.LFEO(1,:) - actBody.LTIO(1,:));
% d_rtibia = norm(actBody.RFEO(1,:) - actBody.RTIO(1,:));
% 
% updateFigureContents('Constraint Check');
% pelib.viz.plotLowerBodySegmentLengthError(estBody, d_pelvis, ...
%     d_lfemur, d_rfemur, d_ltibia, d_rtibia);
% hold on;
% cstrStateL = sum(estState2.cstrStateU(:,1:3), 2) > 0;
% cstrStateLV = zeros(estBody.nSamples, 1);
% cstrStateLV(~cstrStateL) = nan;
% cstrStateR = sum(estState2.cstrStateU(:,4:6), 2) > 0;
% cstrStateRV = zeros(estBody.nSamples, 1);
% cstrStateRV(~cstrStateR) = nan;
% t = 1:estBody.nSamples;
% scatter(t, cstrStateLV, '<c'); scatter(t, cstrStateRV, '>c');
% legend('Hips', 'LFemur', 'RFemur', 'LTibia', 'RTibia', 'LCstr', 'RCstr');

%% Dynamic plots (animation)
az = 0; el = 180;
updateFigureContents('Animation'); 
% tmpBody1 = xsensBody;
tmpBody1 = viconBodyRel;
tmpBody2 = estBodyRel;
tmpBody1Limits = [tmpBody1.xlim() tmpBody1.ylim() tmpBody1.zlim()];
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

% Animation
updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyRel = estBody.changeRefFrame('MIDPEL');
actBodyRel = actBody.changeRefFrame('MIDPEL');
estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
actBody2 = actBodyRel.toWorldFrame(actBody.MIDPEL+[1 0 0], actBody.qRPV);
estBodyLimits = [estBody2.xlim()+[-1 +1] estBody2.ylim()+[-1 +1] estBody2.zlim()];
i = 1; az = 0; el = 180;
while i <= estBody.nSamples
    [az, el] = view;
    clf; grid;
    xlim(estBodyLimits(1:2)); 
    ylim(estBodyLimits(3:4)); 
    zlim(estBodyLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(az, el);
    pelib.viz.plotLowerBody(estBody2, i, true, false);
    pelib.viz.plotLowerBody(actBody2, i, true, false);
    i = i+10;
    pause(1/1000);
end

% Act and Est Body Shadowing
updateFigureContents('Act and Est Shadowing');
i = 1;
actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()];
az = 0; el = 0;
while i <= actBody.nSamples
    [az, el] = view;
    clf; grid;
    xlim(actBodyLimits(1:2)); 
    ylim(actBodyLimits(3:4)); 
    zlim(actBodyLimits(5:6));
    xlabel('x'); ylabel('y'); zlabel('z');
    view(az, el);
    pelib.viz.plotLowerBody(estBody, i, true, false);
    pelib.viz.plotLowerBody(actBody, i, true, false);
    i = i+10;
    pause(1/1000);
end

% Act and Est Rel Body Shadowing
updateFigureContents('Act and Est Rel Shadowing');
actBodyRel = actBody.changeRefFrame('MIDPEL');
estBodyRel = estBody.changeRefFrame('MIDPEL');
actBodyRelLimits = [actBodyRel.xlim() actBodyRel.ylim() actBodyRel.zlim()];
i = 1;
az = 0; el = 0;
while i <= actBody.nSamples
    [az, el] = view;
    clf; grid;
    xlim(actBodyRelLimits(1:2)); 
    ylim(actBodyRelLimits(3:4)); 
    zlim(actBodyRelLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(az, el);
    pelib.viz.plotLowerBody(estBodyRel, i);
    pelib.viz.plotLowerBody(actBodyRel, i);
    i = i+10;
    pause(1/1000);
end

%     updateFigureContents('Animation Freeze');
%     grid; view(0, 90); hold on;
%     for i=idx0(1):30:idx0(end)
%         pelib.viz.plotLowerBody(estBody, i);
%         pelib.viz.plotLowerBody(actBody, i);
%     end
% %     
% %     updateFigureContents('GFR Acc Diff');
% %     diff_gfr_acc_MP = gfr_acc_MP(1:end-1,:) - gfr_acc_MP_act;
% %     diff_gfr_acc_LA = gfr_acc_LA(1:end-1,:) - gfr_acc_LA_act;
% %     diff_gfr_acc_RA = gfr_acc_RA(1:end-1,:) - gfr_acc_RA_act;
% %     pelib.viz.plotXYZ(diff_gfr_acc_MP, diff_gfr_acc_LA, diff_gfr_acc_RA);
%     
%     

%     
%         
updateFigureContents('Animation');
actBodyRelLimits = [actBodyRel.xlim() actBodyRel.ylim() actBodyRel.zlim()];
i = 1;
while i <= actBody.nSamples
    clf; grid;
    xlim(actBodyRelLimits(1:2)); 
    ylim(actBodyRelLimits(3:4)); 
    zlim(actBodyRelLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(20, 30);
    pelib.viz.plotLowerBody(actBodyRel, i, true, false);
    i = i+5;
    pause(1/1000);
end

updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyLimitsRel = [estBodyRel.xlim() estBodyRel.ylim() estBodyRel.zlim()];
i=10;
while i <= estBodyRel.nSamples
    clf; grid;
    xlim(estBodyLimitsRel(1:2)); 
    ylim(estBodyLimitsRel(3:4)); 
    zlim(estBodyLimitsRel(5:6));  
    view(40, 30);
    pelib.viz.plotLowerBody(estBodyRel, i, true, false);
    i = i+20;
    pause(1/1000);
end
