%% Run experiment for all NeuRA data
% To display the mean of results per label
% row = startsWith(results.name, 'tcd-s4') | startsWith(results.name, 'tcd-s5')
% fresults = results(row, :)
% varfun(@mean, results, 'InputVariables', @isnumeric, 'GroupingVariables', 'label')
% writetable(rTable, '
dir = 'neura';
expDir = sprintf('%s/explore', dir);

dataList = { ...
    struct('subj', 'S00', 'act', 'Trial012'), ...
    struct('subj', 'S00', 'act', 'Trial014'), ...
    struct('subj', 'S00', 'act', 'Trial017'), ...
};

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C49', 'R_LowLeg', '00B40C4A', ... 
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
% 'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...

setups = {
    struct('est', 'ekfv3', 'accData', 's', 'oriData', 's', 'accDataNoise', 0, ...
           'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
    struct('est', 'ekfv3', 'accData', 'v', 'oriData', 'v', 'accDataNoise', 0, ...
           'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
};

for mI = [2]
    for cI = [222]
        setups{end+1} = struct('est', 'ekfv3', ...
           'accData', 's', 'oriData', 's', 'accDataNoise', 0, ...
           'applyMeas', mI, 'applyCstr', cI, 'P', 0.5);
    end
end


for i = 1:length(setups)
    setups{i}.label = getLabel(setups{i});
end
           
dataN = length(dataList);
results = table();

for i = 1:dataN
    n = dataList{i};
    
    name = sprintf('%s-%s-%s', 'neura', n.subj, n.act);
    dataPath = sprintf('%s/mat/%s.mat', dir, name);
    if exist(dataPath, 'file')
        load(dataPath, 'data');
    else
        data = struct('name', name, ...
            'fnameV', sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act), ...
            'fnameS', sprintf('%s/imu/%s-%s', dir, n.subj, n.act), ...
            'fnameX', sprintf('%s/xsens/%s-%s.mat', dir, n.subj, n.act) );
        
        data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);
        data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
        % data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
        data.dataX = false;
        save(dataPath, 'data');
    end
    
    display(sprintf("Data %3d/%3d: %s", i, dataN, data.name));
    r = runNeuRAExperiment(data.dataV, data.dataS, data.dataX, ...
                           data.name, setups, expDir);
    results = [results; struct2table(r)];
end

% Append new results
dataPath = sprintf("%s/results.mat", expDir);
if exist(dataPath, 'file')
    newResults = results;
    load(dataPath);
    [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
    results(ia,:) = [];
    results = [results; newResults];
end
save(sprintf("%s/results.mat", expDir), 'results')

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
    label = sprintf('N%s%s+M%02d+C%03d', aD, setup.oriData, ...
        setup.applyMeas, setup.applyCstr);
end