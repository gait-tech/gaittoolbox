make_it_tight = true;
% subplot = @(m,n,p) explore.subtightplot (m, n, p, [0.025 0.05], [0.125 0.02], [0.11 0.02]);
subplot = @(m,n,p) explore.subtightplot (m, n, p, [0.025 0.05], [0.125 0.02], [0.05 0.01]);
if ~make_it_tight,  clear subplot;  end

addpath('etclib');

setup = struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+Aw__sOw__sIw__v+Sav03+M76+C355");
ns = extractBetween(setup.algo, 1, 3);

load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, setup.file));
load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, setup.file, setup.algo));

if cs.initSrc == 'w__v'
    csActBody = W__viconBody;
elseif cs.initSrc == 'v__v'
    csActBody = V__viconBody;
else
    csActBody = W__xsensBody;
end
csActBodyRel = csActBody.changeRefFrame('MIDPEL');    
estBodyRel = estBody.changeRefFrame('MIDPEL');
idxS = find(allIdx.w__x==allIdx.w__v(1));
idxE = find(allIdx.w__x==allIdx.w__v(end));
xsensBodyRel = W__xsensBody.getSubset(idxS:idxE).changeRefFrame('MIDPEL');

estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
xsensBodyRel2 = xsensBodyRel.toWorldFrame(csActBody.MIDPEL, xsensBodyRel.qRPV);

dPos = estBody2.calcDPos(csActBody2);
dOri = estBody2.calcDOri(csActBody2);
    
estKA = rad2deg(estBody2.calcJointAnglesLKnee());
actKA = rad2deg(csActBody2.calcJointAnglesLKnee());
xsnKA = rad2deg(xsensBodyRel2.calcJointAnglesLKnee());
estHA = rad2deg(estBody2.calcJointAnglesLHip());
actHA = rad2deg(csActBody2.calcJointAnglesLHip());
xsnHA = rad2deg(xsensBodyRel2.calcJointAnglesLHip());
t = (1:size(estKA,1))/estBody.fs;

lw = 4;
fontsize = 24;
spH = 4;
spW = 1;

plotXsens = false;

imageSize = get(0, 'Screensize');
% imageSize(3) = imageSize(3)/2;
F = updateFigureContents('sample knee and hip joint angles');
subplot(spH,spW,1);
p1 = plot(t, estKA(:,2), '-k', t, actKA(:,2), ':k', 'LineWidth', lw); grid;
if plotXsens
    hold on;
    p1 = plot(t, xsnKA(:,2), '-.k', 'LineWidth', lw);
end
ylabel('knee Y ({\circ})'); 
xticks(0:10); xticklabels([]); xlim([t(1),t(end)])
subplot(spH,spW,2);
p2 = plot(t, estHA(:,2), '-k', t, actHA(:,2), ':k', 'LineWidth', lw); grid;
if plotXsens
    hold on;
    p2 = plot(t, xsnHA(:,2), '-.k', 'LineWidth', lw);
end
ylabel('hip Y ({\circ})'); 
xticks(0:10); xticklabels([]); xlim([t(1),t(end)])
subplot(spH,spW,3);
p3 = plot(t, estHA(:,1), '-k', t, actHA(:,1), ':k', 'LineWidth', lw); grid;
if plotXsens
    hold on;
    p3 = plot(t, xsnHA(:,1), '-.k', 'LineWidth', lw);
end
ylabel('hip X ({\circ})');  
xticks(0:10); xticklabels([]); xlim([t(1),t(end)])
subplot(spH,spW,4);
p4 = plot(t, estHA(:,3), '-k', t, actHA(:,3), ':k', 'LineWidth', lw); grid;
if plotXsens
    hold on;
    p4 = plot(t, xsnHA(:,3), '-.k', 'LineWidth', lw);
end
xlabel('Time (s)');  ylabel('hip Z ({\circ})');
xticks(0:10);
xlim([t(1),t(end)])
if plotXsens
    lgh = legend('CKF', 'Vicon', 'Xsens', 'Orientation', 'horizontal');
else
    lgh = legend('CKF', 'Vicon', 'Orientation', 'horizontal');
end
% set(lgh, 'position', [0.17    0.025    0.13    0.0496]);
set(lgh, 'position', [0.05    0.025    0.13    0.0496]);
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
set(F, 'Position', imageSize);
saveas(F, 'results-kneehip-angle-viconsample4.png', 'png');


%set(findall(gcf,'-property','XMinorTick'),'XMinorTick', 'on')
%set(findall(gcf,'-property','YMinorTick'),'YMinorTick', 'on')

% fontsize = 20;
% spH = 3;
% spW = 1;
% clf; updateFigureContents('knee and hip joint angles + dpos and dori');
% subplot(spH,spW,1);
% p1 = plot(t, estKA(:,2), '-b', t, actKA(:,2), ':r', 'LineWidth', lw); grid;
% xlabel('Time (s)'); ylabel('knee Y ({\circ})'); 
% xticks(0:10); xticklabels([]);
% lgh = legend('CKF', 'Vicon', 'Orientation', 'horizontal');
% subplot(spH,spW,2);
% p5 = plot(t, dPos, '-b', 'LineWidth', lw); grid;
% xlabel('Time (s)');  ylabel('$e_{pos}$ (m)', 'Interpreter', 'latex');
% xticks(0:10);
% subplot(spH,spW,3);
% p6 = plot(t, dOri, '-b', 'LineWidth', lw); grid;
% xlabel('Time (s)');  ylabel('$e_{ori}$ ($^\circ$)', 'Interpreter', 'latex');
% xticks(0:10);
% % set(lgh, 'position', [0.17    0.00    0.13    0.0496]);
% set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
% set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')

% imageSize = get(0, 'Screensize');
% fontsize = 32;
% 
% subplot = @(m,n,p) explore.subtightplot (m, n, p, [0.06 0.02], [0.12 0.02], [0.06 0.01]);
% 
% %epos and eori
% T = readtable("+explore/tech1-eposeori.csv");
% F = updateFigureContents('epos eori');
% idx=1:5; subplot(2,1,1); 
% bar(categorical(T.x(idx), T.x(idx)), [T.Vicon(idx) T.Xsens(idx) T.CKF(idx)]);
% xlabel('Types of Motion'); 
% ylabel('$e_{pos}$ (cm)', 'Interpreter', 'latex'); 
% tmpy = yticklabels; tmpy(1) = {'  0'}; yticklabels(tmpy);
% grid on;
% 
% idx=6:10; subplot(2,1,2); 
% bar(categorical(T.x(idx), T.x(idx)), [T.Vicon(idx) T.Xsens(idx) T.CKF(idx)]);
% xlabel('Types of Motion'); 
% ylabel('$e_{ori}$ ($^\circ$)', 'Interpreter', 'latex'); 
% lgh = legend('\it{Vicon input}', '\it{Xsens}', '\it{CKF}', ...
%        'Location', 'north', 'Orientation', 'horizontal');
% grid on;
% set(lgh, 'position', [0.13    0.01    0.13    0.0496]);
% set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
% set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
% set(findall(gcf,'-property','LineWidth'),'LineWidth', 2)
% set(F, 'Position', imageSize);
% % applyhatch(gcf,'\-x.');
% % h,patterns,CvBW,Hinvert,colorlist,dpi,hatchsc,linewidth
% % [im_hatch,colorlist] = applyhatch_pluscolor(gcf,'wk\',0,1,[],100,3,12);
% saveas(F, 'results-dposdorimean.png', 'png');
% % saveas(gcf, 'results-dposdorimean-bw.png', 'png');
% 
% %rmse and cc
% T = readtable("+explore/tech1-jointangles.csv");
% F = updateFigureContents('joint rmse cc');
% idx=1:5; subplot(2,1,1); 
% bar(categorical(T.x(idx), T.x(idx)), [T.kneeY(idx) T.hipY(idx) T.hipX(idx) T.hipZ(idx)]);
% xlabel('Types of Motion'); ylabel('Joint Angle RMSE ({\circ})'); 
% grid on;
% idx=6:10; subplot(2,1,2); 
% bar(categorical(T.x(idx), T.x(idx)), [T.kneeY(idx) T.hipY(idx) T.hipX(idx) T.hipZ(idx)]);
% xlabel('Types of Motion'); ylabel('Joint Angle CC'); ylim([0 1]);
% lgh = legend('knee Y', 'hip Y', 'hip X', 'hip Z', ...
%        'Location', 'northwest', 'Orientation', 'horizontal');
% grid on;
% set(lgh, 'position', [0.15    0.01    0.13    0.0496]);
% set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
% set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
% set(findall(gcf,'-property','LineWidth'),'LineWidth', 2)
% % [im_hatch,colorlist] = applyhatch_pluscolor(gcf,'k\wxk',0,0,[],100,5,3);
% set(F, 'Position', imageSize);
% saveas(F, 'results-kneehiprmsecc.png', 'png');
% % saveas(gcf, 'results-kneehiprmsecc-bw.png', 'png');