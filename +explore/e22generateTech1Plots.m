make_it_tight = true;
subplot = @(m,n,p) explore.subtightplot (m, n, p, [0.02 0.05], [0.105 0.02], [0.11 0.02]);
if ~make_it_tight,  clear subplot;  end

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
estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);

dPos = estBody2.calcDPos(csActBody2);
dOri = estBody2.calcDOri(csActBody2);
    
estKA = rad2deg(estBody2.calcJointAnglesLKnee());
actKA = rad2deg(csActBody2.calcJointAnglesLKnee());
estHA = rad2deg(estBody2.calcJointAnglesLHip());
actHA = rad2deg(csActBody2.calcJointAnglesLHip());
t = (1:size(estKA,1))/estBody.fs;

lw = 4;
fontsize = 24;
spH = 4;
spW = 1;

clf; updateFigureContents('sample knee and hip joint angles');
subplot(spH,spW,1);
p1 = plot(t, estKA(:,2), '-b', t, actKA(:,2), ':r', 'LineWidth', lw); grid;
xlabel('Time (s)'); ylabel('knee Y ({\circ})'); 
xticks(0:10); xticklabels([]);
subplot(spH,spW,2);
p2 = plot(t, estHA(:,2), '-b', t, actHA(:,2), ':r', 'LineWidth', lw); grid;
xlabel('Time (s)'); ylabel('hip Y ({\circ})'); 
xticks(0:10); xticklabels([]);
subplot(spH,spW,3);
p3 = plot(t, estHA(:,1), '-b', t, actHA(:,1), ':r', 'LineWidth', lw); grid;
xlabel('Time (s)'); ylabel('hip X ({\circ})');  
xticks(0:10); xticklabels([]);
subplot(spH,spW,4);
p4 = plot(t, estHA(:,3), '-b', t, actHA(:,3), ':r', 'LineWidth', lw); grid;
xlabel('Time (s)');  ylabel('hip Z ({\circ})');
xticks(0:10);
lgh = legend('CKF', 'Vicon', 'Orientation', 'horizontal');
set(lgh, 'position', [0.17    0.00    0.13    0.0496]);
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
%set(findall(gcf,'-property','XMinorTick'),'XMinorTick', 'on')
%set(findall(gcf,'-property','YMinorTick'),'YMinorTick', 'on')

fontsize = 20;
spH = 3;
spW = 1;
clf; updateFigureContents('knee and hip joint angles + dpos and dori');
subplot(spH,spW,1);
p1 = plot(t, estKA(:,2), '-b', t, actKA(:,2), ':r', 'LineWidth', lw); grid;
xlabel('Time (s)'); ylabel('knee Y ({\circ})'); 
xticks(0:10); xticklabels([]);
lgh = legend('CKF', 'Vicon', 'Orientation', 'horizontal');
subplot(spH,spW,2);
p5 = plot(t, dPos, '-b', 'LineWidth', lw); grid;
xlabel('Time (s)');  ylabel('$e_{pos}$ (m)', 'Interpreter', 'latex');
xticks(0:10);
subplot(spH,spW,3);
p6 = plot(t, dOri, '-b', 'LineWidth', lw); grid;
xlabel('Time (s)');  ylabel('$e_{ori}$ ($^\circ$)', 'Interpreter', 'latex');
xticks(0:10);
% set(lgh, 'position', [0.17    0.00    0.13    0.0496]);
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')

imageSize = get(0, 'Screensize');
fontsize = 32;

subplot = @(m,n,p) explore.subtightplot (m, n, p, [0.06 0.02], [0.12 0.02], [0.06 0.01]);

%epos and eori
T = readtable("+explore/tech1-eposeori.csv");
F = updateFigureContents('epos eori'); set(F, 'Position', imageSize);
idx=1:5; subplot(2,1,1); 
bar(categorical(T.x(idx), T.x(idx)), [T.Vicon(idx) T.Xsens(idx) T.CKF(idx)])
xlabel('Types of Motion'); 
ylabel('$e_{pos}$ (cm)', 'Interpreter', 'latex'); 
tmpy = yticklabels; tmpy(1) = {'  0'}; yticklabels(tmpy);
grid on;

idx=6:10; subplot(2,1,2); 
bar(categorical(T.x(idx), T.x(idx)), [T.Vicon(idx) T.Xsens(idx) T.CKF(idx)])
xlabel('Types of Motion'); 
ylabel('$e_{ori}$ ($^\circ$)', 'Interpreter', 'latex'); 
lgh = legend('\it{Vicon input}', '\it{Xsens}', '\it{CKF}', ...
       'Location', 'north', 'Orientation', 'horizontal');
grid on;
set(lgh, 'position', [0.13    0.01    0.13    0.0496]);
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
saveas(F, 'results-dposdorimean.png', 'png');

%rmse and cc
T = readtable("+explore/tech1-jointangles.csv");
F = updateFigureContents('joint rmse cc'); set(F, 'Position', imageSize);
idx=1:5; subplot(2,1,1); 
bar(categorical(T.x(idx), T.x(idx)), [T.kneeY(idx) T.hipY(idx) T.hipX(idx) T.hipZ(idx)])
xlabel('Types of Motion'); ylabel('Joint Angle RMSE ({\circ})'); 
grid on;
idx=6:10; subplot(2,1,2); 
bar(categorical(T.x(idx), T.x(idx)), [T.kneeY(idx) T.hipY(idx) T.hipX(idx) T.hipZ(idx)])
xlabel('Types of Motion'); ylabel('Joint Angle CC'); ylim([0 1]);
lgh = legend('knee Y', 'hip Y', 'hip X', 'hip Z', ...
       'Location', 'northwest', 'Orientation', 'horizontal');
grid on;
set(lgh, 'position', [0.15    0.01    0.13    0.0496]);
set(findall(gcf,'-property','FontSize'),'FontSize', fontsize)
set(findall(gcf,'-property','FontName'),'FontName', 'Times New Roman')
saveas(F, 'results-kneehiprmsecc.png', 'png');