dir = 'neura-sparse01';
dataList = { ...
    struct('subj', 'S01', 'act', 'Trial-Walk-2'), ...
};

subject = {'S01'};
target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'Jog-1', 'Jog-2', ...
    'JumpingJacks-1', 'JumpingJacks-2', 'FigureofEight-1', 'FigureofEight-2', ...
    'Zigzag-1', 'Zigzag-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
    'HighKneeJog-1', 'HighKneeJog-2', 'Fivemin-1', 'Fivemin-2'};
for i=1:length(target)
    for s=1:length(subject)    
        dataList{end+1} = struct('subj', subject{s}, 'act', sprintf("Trial-%s", target{i}));
    end
end
 
subject = {'S02', 'S03', 'S04', 'S05', 'S06', 'S07', 'S08', 'S09', 'S10'};
target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
    'JumpingJacks-1', 'JumpingJacks-2', 'FigureofEight-1', 'FigureofEight-2', ...
    'Zigzag-1', 'Zigzag-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
    'HighKneeJog-1', 'HighKneeJog-2', 'Fivemin-1', 'Fivemin-2'};
for i=1:length(target)
    for s=1:length(subject)    
        dataList{end+1} = struct('subj', subject{s}, 'act', sprintf('Trial-%s', target{i}));
    end
end

T = table();
dataN = length(dataList);
for i = 1:dataN
    n = dataList{i};
    
    name = sprintf("%s-%s-%s", 'neura', n.subj, n.act);
    fnameV = sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act);
    dataV = mocapdb.ViconBody.loadViconMat(fnameV);
    n.startFrame = min(findMovementFromStartFrame(dataV.LTIO), ...
                       findMovementFromStartFrame(dataV.RTIO));
    n.endFrame = max(findMovementFromEndFrame(dataV.LTIO), ...
                     findMovementFromEndFrame(dataV.RTIO));
    
    dataList{i} = n;
    
    LTIO = dataV.LTIO; RTIO = dataV.RTIO;
    minP = min([min(LTIO); min(RTIO)]);
    maxP = max([max(LTIO); max(RTIO)]);
    clf; pelib.viz.plotXYZ(1, LTIO, RTIO); hold on;
    for j=1:3
        subplot(3, 1, j);
        line([n.startFrame n.startFrame], [minP(j) maxP(j)]);
        line([n.endFrame n.endFrame], [minP(j) maxP(j)]);
    end
    
    T = [T; struct2table(n)];
    saveas(gcf, sprintf('explore_output/StartEndFrame/%s.png', name))
end

writetable(T, sprintf('%s/data-list.csv', dir));

function pos = findMovementFromStartFrame(point)
    d = vecnorm(point - point(1, :), 2, 2);
    pos = find(d > 10);
    if length(pos) > 0
        pos = max(pos(1)-10, 0);
    else
        pos = -1;
    end
end
function pos = findMovementFromEndFrame(point)
    d = vecnorm(point - point(end, :), 2, 2);
    pos = find(d > 10);
    if length(pos) > 0
        pos = max(pos(end)+10, 0);
    else
        pos = -1;
    end
end