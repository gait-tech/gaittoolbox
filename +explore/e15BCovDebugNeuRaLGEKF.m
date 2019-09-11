% motion list
ds = struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P101+M125+C007");

dataSfname = sprintf('neura-sparse01/imu/%s', ds.file);
ns = extractBetween(ds.algo, 1, 3);

load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, ds.file));
load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, ds.file, ds.algo));

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
sensors = {};
[estStateDebug, estPDebug] = experiment.buildStateDebug(estState, estState2, cs.est);
[estBodyDebug, sensorsDebug] = experiment.buildgrBodyDebug(estBody, sensors, estState2, cs.est);
[viconBodyDebug, ~] = experiment.buildgrBodyDebug(vb, {}, false, 'vicon');
   
%% Animation
% az = 0; el = 180;
updateFigureContents('Animation'); 
tmpBody1 = estBodyDebug;
tmpBody2 = false;
tmpBody1Limits = [tmpBody1.xlim() tmpBody1.ylim() tmpBody1.zlim()];
samplePlotN = 1000;
% for i=3:3*5:estBodyDebug.nSamples
for i=1
    [az, el] = view;
    clf; grid; axis equal;
    xlabel('x'); ylabel('y'); zlabel('z');
    xlim(tmpBody1Limits(1:2)); 
    ylim(tmpBody1Limits(3:4)); 
    zlim(tmpBody1Limits(5:6));  
    view(az, el);
    pelib.viz.plotLowerBody(tmpBody1, i, true, false);
    hold on;
    
    ps = explore.e15sample(estStateDebug, estPDebug, samplePlotN, i, 1.0/estBody.fs, cs.est);
    scatter3(ps.PV(:,1), ps.PV(:,2), ps.PV(:,3), 'g.');
    scatter3(ps.LS(:,1), ps.LS(:,2), ps.LS(:,3), 'b.');
    scatter3(ps.RS(:,1), ps.RS(:,2), ps.RS(:,3), 'r.');
%     pelib.viz.plotMeanCovSamples(tmpBody1.RTIO(i, :), ...
%                 estPRAPosDebug(:, :, i), samplePlotN, ...
%                 sprintf('%s.', tmpBody1.xyzColor{1}));
%     pelib.viz.plotMeanCovSamples(tmpBody1.MIDPEL(i, :), ...
%                 estPMPPosDebug(:, :, i), samplePlotN, ...
%                 sprintf('%s.', tmpBody1.xyzColor{2}));
%     pelib.viz.plotMeanCovSamples(tmpBody1.LTIO(i, :), ...
%                 estPLAPosDebug(:, :, i), samplePlotN, ...
%                 sprintf('%s.', tmpBody1.xyzColor{3}));
    if ~(tmpBody2==false)
        pelib.viz.plotLowerBody(tmpBody2, i, true, false);
    end
    
    pause(1/1000);
end