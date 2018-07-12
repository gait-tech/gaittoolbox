%% Run experiment for all NeuRA data
% To display the mean of results per label
% row = startsWith(results.name, 'tcd-s4') | startsWith(results.name, 'tcd-s5')
% fresults = results(row, :)
% varfun(@mean, results, 'InputVariables', @isnumeric, 'GroupingVariables', 'label')
% writetable(rTable, '
dir = 'neura';
expDir = sprintf('%s/explore', dir);

dataList = { ...
    struct('subj', 'S00-2', 'act', 'Trial-016'), ...
    struct('subj', 'S00-2', 'act', 'Trial-018'), ...
%     struct('subj', 'S00-2', 'act', 'Trial-020'), ...
%     struct('subj', 'S00-2', 'act', 'Trial-023'), ...
%     struct('subj', 'S00-2', 'act', 'Trial-024'), ...
%     struct('subj', 'S00-2', 'act', 'Trial-025'), ...
%     struct('subj', 'S00-2', 'act', 'Trial-026'), ...
};

options = struct('Pelvis', 'Pelvis', ...
    'L_UpLeg', 'LeftUpperLeg', 'R_UpLeg', 'RightUpperLeg', ...
    'L_LowLeg', 'prop', 'R_LowLeg', 'prop_1', ...
    'L_Foot', 'LeftFoot', 'R_Foot', 'RightFoot');

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
            'fnameX', sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act), ...
            'fnameS', sprintf('%s/imu/%s-%s.mvnx', dir, n.subj, n.act));
        
        data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);
        data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
        data.dataS = mocapdb.XsensBody.loadMVNX(data.fnameS, options);
        save(dataPath, 'data');
    end
    
    display(sprintf("Data %3d/%3d: %s", i, dataN, data.name));
    r = runNeuRAExperiment(data.dataV, data.dataX, data.dataS, ...
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