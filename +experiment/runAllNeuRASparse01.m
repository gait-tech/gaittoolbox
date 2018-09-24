% ======================================================================
%> Run experiment for all NeuRA data
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore', dir);

dataList = { ...
%     struct('subj', 'S01', 'act', 'Trial-Walk-1')
    struct('subj', 'S02', 'act', 'Trial-Walk-1')
    struct('subj', 'S03', 'act', 'Trial-Walk-1')
};
% subject = {'S01', 'S02', 'S03', 'S04'};
% target = {'Static-1', 'Walk-1', 'Walk-2', 'TUG-1', 'TUG-2', 'Jog-1', 'Jog-2', ...
%     'JumpingJacks-1', 'JumpingJacks-2', 'FigureofEight-1', 'FigureofEight-2', ...
%     'Zigzag-1', 'Zigzag-2', 'SpeedSkater-1', 'SpeedSkater-2', ...
%     'HighKneeJog-1', 'HighKneeJog-2', 'Fivemin-1', 'Fivemin-2'};
% for i=1:length(target)
%     for s=1:length(subject)    
%         dataList{end+1} = struct('subj', subject{s}, 'act', sprintf('Trial-%s', target{i}));
%     end
% end

% options = struct('Pelvis', '00B40B91', ...
%     'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%     'L_LowLeg', '00B40C49', 'R_LowLeg', '00B40C4A', ...
%     'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
%     'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
% options = struct('Pelvis', 'Pelvis', ...
%     'L_UpLeg', 'LeftUpperLeg', 'R_UpLeg', 'RightUpperLeg', ...
%     'L_LowLeg', 'prop', 'R_LowLeg', 'prop_1', ...
%     'L_Foot', 'LeftFoot', 'R_Foot', 'RightFoot');

setups = {
%     struct('est', 'ekfv3', 'accData', 's', 'oriData', 's', 'initSrc', 'v', ...
%            'accDataNoise', 0, 'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
%     struct('est', 'ekfv3', 'accData', 'v', 'oriData', 'v', 'accDataNoise', 0, ...
%            'applyMeas', 21, 'applyCstr', 0, 'P', 0.5), ...
};

% for mI = [0]
%     for cI = [0 1:8 21:23 51:54 71:77 121:122 124:125 131:132 134:135 ...
%               141:144 151:154 201:208 221:223 271:278]
%         setups{end+1} = struct('est', 'ekfv3', ...
%            'accData', 'v', 'oriData', 'v', 'accDataNoise', 0, ...
%            'initSrc', 'v', 'applyMeas', mI, 'applyCstr', cI, 'P', 0.5);
%     end
% end

for iI = {'w__v'} % ['v' 'x']
    for mI = [2]
%         for cI = [0 1:8 21:23 51:54 71:78 121:122 124:125 131:132 134:135 ...
%               141:144 151:154 201:208 221:223 271:278]
%         for cI = [0 201:208 221:223 132 135 151:154 ]
%         for cI = [202 203 135 152]
        for cI = [201 202 203 152 172]
%         for cI = [0 1:8 21:23 122 125 141:144 ...
%                     201:208 221:223 132 135 151:154 ]
            for aI = {'w__s'}
                for oI = {'w__s'}
                    setups{end+1} = struct('est', 'ekfv3', 'frame', 'world', ...
                       'accData', aI, 'oriData', oI, 'accDataNoise', 0, ...
                       'initSrc', iI, 'applyMeas', mI, 'applyCstr', cI, ...
                       'P', 0.5);
                end
            end
        end
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
            'fnameS', sprintf('%s/imu/%s-%s', dir, n.subj, n.act));
        
        data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);
        if exist(data.fnameX, 'file')
            data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
        else
            data.dataX = [];
        end
%         data.dataS = mocapdb.XsensBody.loadMVNX(data.fnameS, options);
        data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
        data.dataS.fs = 100;
        data.calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
             sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
             sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );
        data.calibW2V = mocapdb.XsensBody('nSamples', 0);
        save(dataPath, 'data');
    end
    
    display(sprintf("Data %3d/%3d: %s", i, dataN, data.name));
    r = experiment.runNeuRASparse01Experiment(data.dataS, ...
            data.dataV, data.calibV2W, data.calibW2V, ...
            data.dataX, ...
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
    label = sprintf('NS1+A%sO%sI%s+M%02d+C%03d', aD, setup.oriData, setup.initSrc, ...
        setup.applyMeas, setup.applyCstr);
end