%% Test vicon load
dataV = mocapdb.ViconBody.loadViconMat('+unittest/viconBodyData/S03-Trial-Walk-1.mat');
r = [9.654742,-0.718345,2.264074,-0.047825,0.060515,-0.021490,0.994141,-0.391846,0.300049,0.392628,0.516707,-0.339639,0.680810];
assert(all(dataV.LFEP(28, :) == [29.8426000000000,1338,804.732000000000]));
assert(all(dataV.RFEO(14, :) == [-145.438000000000,1380.89000000000,456.712000000000]));
assert(all(dataV.LTIO(4, :) == [25.3427000000000,1378.59000000000,88.7196000000000]));

%% Test subset data
dataV = mocapdb.ViconBody.loadViconMat('+unittest/viconBodyData/S03-Trial-Walk-1.mat');
idx = 1:3:300;
dataVSubset = dataV.getSubset(idx);

segList = [dataV.posList dataV.oriList];
for i=1:length(segList)
    n = segList{i};
    if sum(size(dataV.(n))) ~= 0
        assert(all(all(dataVSubset.(n) == dataV.(n)(idx, :))));
    end
end
assert(dataVSubset.nSamples == length(idx));

%% Test change pos unit
dataV = mocapdb.ViconBody.loadViconMat('+unittest/viconBodyData/S03-Trial-Walk-1.mat');
dataV2 = dataV.changePosUnit('m');

for i=1:length(dataV.posList)
    n = dataV.posList{i};
    if sum(size(dataV.(n))) ~= 0
        assert(all(all((dataV2.(n)-dataV.(n)/1000)<1e-10)));
    end
end
assert(dataV2.posUnit == 'm');