% Show the absolute and relative position
fs=100;
updateFigureContents('Hips');
clf; grlib.viz.plotStateComparison(estState2, actState, 1, fs);

updateFigureContents('LAnkle');
clf; grlib.viz.plotStateComparison(estState2, actState, 11, fs);

updateFigureContents('RAnkle');
clf; grlib.viz.plotStateComparison(estState2, actState, 21, fs);

estStateRel2 = grlib.est.changeStateRefFrame(estState2);
actStateRel = grlib.est.changeStateRefFrame(actState);
updateFigureContents('LAnkleRel');
clf; grlib.viz.plotStateComparison(estStateRel2, actStateRel, 11, 60);

updateFigureContents('RAnkleRel');
clf; grlib.viz.plotStateComparison(estStateRel2, actStateRel, 21, 60);

updateFigureContents('Hips2');
clf; grlib.viz.plotStateComparison(estState2, actState, 1, 60);

updateFigureContents('LAnkle2');
clf; grlib.viz.plotStateComparison(estState2, actState, 11, 60);

updateFigureContents('RAnkle2');
clf; grlib.viz.plotStateComparison(estState2, actState, 21, 60);

estStateRel2 = grlib.est.changeStateRefFrame(estState2);
actStateRel = grlib.est.changeStateRefFrame(actState);
updateFigureContents('LAnkleRel2');
clf; grlib.viz.plotStateComparison(estStateRel2, actStateRel, 11, 60);

updateFigureContents('RAnkleRel2');
clf; grlib.viz.plotStateComparison(estStateRel2, actStateRel, 21, 60);

% P exploration
load('C:\Users\z5151460\workspace\gaitrecon\experiments\tcd-s1-acting1-Dxxx+ZUPT+C001.mat')
[nSamples, nState] = size(estState);

P_pri1 = nan(nSamples, nState);
P_pos1 = nan(nSamples, nState);
P_con1 = nan(nSamples, nState);

for i=1:nSamples
    P_pri1(i,:) = diag(estState2.predP(:,:,i));
    P_pos1(i,:) = diag(estState2.zuptP(:,:,i));
    P_con1(i,:) = diag(estState2.cstrP(:,:,i));
end

load('C:\Users\z5151460\workspace\gaitrecon\experiments\tcd-s1-acting1-Dxxx+ZUPT+C002.mat')
[nSamples, nState] = size(estState);

P_pri2 = nan(nSamples, nState);
P_pos2 = nan(nSamples, nState);
P_con2 = nan(nSamples, nState);

for i=1:nSamples
    P_pri2(i,:) = diag(estState2.predP(:,:,i));
    P_pos2(i,:) = diag(estState2.zuptP(:,:,i));
    P_con2(i,:) = diag(estState2.cstrP(:,:,i));
end

% Static Plots
grlib.viz.plotPosition({actBody}, {'LTIO', 'RTIO'});

updateFigureContents('Position1');
actBodyRel = actBody.changeRefFrame('MIDPEL');
estBodyRel = estBody.changeRefFrame('MIDPEL');
grlib.viz.plotPosition({estBodyRel, actBodyRel}, {'LTIO', 'RTIO'});

updateFigureContents('Position2');
actBodyRel = actBody.changeRefFrame('MIDPEL');
estBodyRel = estBody.changeRefFrame('MIDPEL');
grlib.viz.plotPosition({estBodyRel, actBodyRel}, {'LTIO', 'RTIO'});

% Constraint Info
d_pelvis = norm(actBody.LFEP(1,:) - actBody.RFEP(1,:));
d_lfemur = norm(actBody.LFEP(1,:) - actBody.LFEO(1,:));
d_rfemur = norm(actBody.RFEP(1,:) - actBody.RFEO(1,:));
d_ltibia = norm(actBody.LFEO(1,:) - actBody.LTIO(1,:));
d_rtibia = norm(actBody.RFEO(1,:) - actBody.RTIO(1,:));

updateFigureContents('Constraint Check');
grlib.viz.plotLowerBodySegmentLengthError(estBody, d_pelvis, ...
    d_lfemur, d_rfemur, d_ltibia, d_rtibia);
hold on;
cstrStateL = sum(estState2.cstrStateU(:,1:3), 2) > 0;
cstrStateLV = zeros(estBody.nSamples, 1);
cstrStateLV(~cstrStateL) = nan;
cstrStateR = sum(estState2.cstrStateU(:,4:6), 2) > 0;
cstrStateRV = zeros(estBody.nSamples, 1);
cstrStateRV(~cstrStateR) = nan;
t = 1:estBody.nSamples;
scatter(t, cstrStateLV, '<c'); scatter(t, cstrStateRV, '>c');
legend('Hips', 'LFemur', 'RFemur', 'LTibia', 'RTibia', 'LCstr', 'RCstr');

% Animation
actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()];
i = 1;
while i <= actBody.nSamples
    clf; grid;
    xlim(actBodyLimits(1:2)); 
    ylim(actBodyLimits(3:4)); 
    zlim(actBodyLimits(5:6));  
    xlabel('x'); ylabel('y'); zlabel('z');
    view(20, 30);
    grlib.viz.plotLowerBody(actBody, i);
    i = i+2;
    pause(1/1000);
end

% Animation
updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyLimits = [estBody.xlim() estBody.ylim() estBody.zlim()];
i = 1; az = 0; el = 180;
while i <= estBody.nSamples
    [az, el] = view;
    clf; grid;
    xlim(estBodyLimits(1:2)); 
    ylim(estBodyLimits(3:4)); 
    zlim(estBodyLimits(5:6));  
    view(az, el);
    grlib.viz.plotLowerBody(estBody, i);
    i = i+1;
    pause(1/1000);
end

% Act and Est Body Shadowing
updateFigureContents('Act and Est Shadowing');
i = 1;
az = 0; el = 0;
while i <= actBody.nSamples
    [az, el] = view;
    clf; grid;
    xlabel('x'); ylabel('y'); zlabel('z');
    view(az, el);
    grlib.viz.plotLowerBody(estBody, i, true);
    grlib.viz.plotLowerBody(actBody, i, true);
    i = i+10;
    pause(1/1000);
end

%     updateFigureContents('Animation Freeze');
%     grid; view(0, 90); hold on;
%     for i=idx0(1):30:idx0(end)
%         grlib.viz.plotLowerBody(estBody, i);
%         grlib.viz.plotLowerBody(actBody, i);
%     end
% %     
% %     updateFigureContents('GFR Acc Diff');
% %     diff_gfr_acc_MP = gfr_acc_MP(1:end-1,:) - gfr_acc_MP_act;
% %     diff_gfr_acc_LA = gfr_acc_LA(1:end-1,:) - gfr_acc_LA_act;
% %     diff_gfr_acc_RA = gfr_acc_RA(1:end-1,:) - gfr_acc_RA_act;
% %     grlib.viz.plotXYZ(diff_gfr_acc_MP, diff_gfr_acc_LA, diff_gfr_acc_RA);
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
    grlib.viz.plotLowerBody(actBodyRel, i);
    i = i+2;
    pause(1/1000);
end

updateFigureContents('Animation');
xlabel('x'); ylabel('y'); zlabel('z');
estBodyLimitsRel = [estBodyRel.xlim() estBodyRel.ylim() estBodyRel.zlim()];
for i=1:estBodyRel.nSamples
    clf; grid;
    xlim(estBodyLimitsRel(1:2)); 
    ylim(estBodyLimitsRel(3:4)); 
    zlim(estBodyLimitsRel(5:6));  
    view(40, 30);
    grlib.viz.plotLowerBody(estBodyRel, i);
    pause(1/1000);
end
