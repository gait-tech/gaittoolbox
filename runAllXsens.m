%% Run experiment to all TCD data
% To display the mean of results per label
% row = startsWith(results.name, 'tcd-s4') | startsWith(results.name, 'tcd-s5')
% fresults = results(row, :)
% varfun(@mean, results, 'InputVariables', @isnumeric, 'GroupingVariables', 'label')
% writetable(rTable, '
dir = 'xsens';
expDir = "./eXsens";

dataList = { ...
    struct('subj', 's1', 'act', 'MT_012000EB-000'), ...
    struct('subj', 's1', 'act', 'MT_012000EB-001'), ...
    struct('subj', 's1', 'act', 'MT_012000EB-002'), ...
    struct('subj', 's1', 'act', 'MT_012000EB-003'), ...
    struct('subj', 's1', 'act', 'MT_012000EB-004'), ...
};

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

setups = {
    struct('est', 'ekfv3', 'accData', 'raw', 'oriData', 'x', ...
           'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
    struct('est', 'ekfv3', 'accData', 'sim', 'oriData', 'x', ...
           'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
};
for mI = [0 1:4]
    for cI = [0 1:8 21:23 51:54 71:78 201:208 221:223 271:278]
        setups{end+1} = struct('est', 'ekfv3', ...
           'accData', 'raw', 'oriData', 'x', ...
           'applyMeas', mI, 'applyCstr', cI, 'P', 0.5);
        setups{end+1} = struct('est', 'ekfv3', ...
           'accData', 'sim', 'oriData', 'x', ...
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
    name = sprintf('%s-%s-%s', 'xsens', n.subj, lower(n.act));
    dataPath = sprintf('%s/mat/%s.mat', dir, name);
    if exist(dataPath, 'file')
        load(dataPath, 'data');
    else
        data = struct('name', name, ...
            'fnameBVH', sprintf('%s/%s.bvh', dir, n.act), ...
            'fnameRaw', sprintf('%s/%s', dir, n.act));
        data.dataBVH = tcdlib.BVHBody.loadXsensBVHFile(data.fnameBVH, "mm");
        data.dataRaw = tcdlib.XsensBody.loadMTExport(data.fnameRaw, options);
        save(dataPath, 'data');
    end
    display(sprintf("Data %3d/%3d: %s", i, dataN, data.name));
    r = runXsensExperiment(data.dataBVH, data.dataRaw, data.name, setups, expDir);
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
    if setup.accData == 'sim'
        if ~isfield(setup, 'accDataNoise') || setup.accDataNoise == 0 
            aD = 's';
        else
            aD = strrep(sprintf('s%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData(1);
    end
    label = sprintf('D%s%s+M%02d+C%03d', aD, setup.oriData, ...
        setup.applyMeas, setup.applyCstr);
end