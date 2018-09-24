%% Recalculate neura error from saved states

dir = 'neura';
expDir = sprintf('%s/explore', dir);

dataList = { ...
};
subject = {'S01', 'S02', 'S03', 'S04'};
target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
    'JumpingJacks-1', 'JumpingJacks-2', };
for i=1:length(target)
    for s=1:length(subject)
        dataList{end+1} = struct('subj', subject{s}, 'act', sprintf('Trial-%s', target{i}));
    end
end

subject = {'S01', 'S02', 'S03'};
target = {'FigureofEight-1', 'FigureofEight-2', 'Zigzag-1', 'Zigzag-2', ...
    'SpeedSkater-1', 'SpeedSkater-2', 'HighKneeJog-1', 'HighKneeJog-2', 'Fivemin-1', 'Fivemin-2'};
for i=1:length(target)
    for s=1:length(subject)
        dataList{end+1} = struct('subj', subject{s}, 'act', sprintf('Trial-%s', target{i}));
    end
end

setups = {};

for iI = ['v']
    for mI = [2]
        for cI = [202 203 135 152 ]
            setups{end+1} = struct('est', 'ekfv3', ...
               'accData', 's', 'oriData', 's', 'accDataNoise', 0, ...
               'initSrc', iI, 'applyMeas', mI, 'applyCstr', cI, ...
               'P', 0.5);
        end
    end
end

for i = 1:size(setups, 2)
    setups{i}.label = getLabel(setups{i});
end

resultsIdx = 1; clear results;
for i = 1:size(dataList, 2)
    n = dataList{i};
    
    name = sprintf('%s-%s-%s', 'neura', n.subj, n.act);
    load(sprintf('%s/%s-debug', expDir, name));
    for j=1:size(setups, 2)
        load(sprintf('%s/%s-%s', expDir, name, setups{j}.label));

        if setups{j}.initSrc == 'v'
            viconBodyRel = viconBody.changeRefFrame('MIDPEL');
            csActBodyRel = viconBodyRel;
        else
            xsensBodyRel = xsensBody.changeRefFrame('MIDPEL');
            csActBodyRel = xsensBodyRel;
        end
        estBodyRel = estBody.changeRefFrame('MIDPEL');
        results0 = estBodyRel.diffRMSE(csActBodyRel);
        
        results0.name = name;
        results0.label = setups{j}.label;
        results(resultsIdx) = results0;
        display(sprintf('%s/%s-%s', expDir, name, setups{j}.label));
        resultsIdx = resultsIdx + 1;
    end
end

results = struct2table(results);

function label = getLabel(setup)
    if setup.accData == 'v'
        if setup.accDataNoise == 0 
            aD = 'v';
        else
            aD = strrep(sprintf('v%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData;
    end
    label = sprintf('N%s%s%s+M%02d+C%03d', aD, setup.oriData, setup.initSrc, ...
        setup.applyMeas, setup.applyCstr);
end