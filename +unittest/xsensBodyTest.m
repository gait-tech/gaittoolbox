%% Test loadMTExport
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataS = mocapdb.XsensBody.loadMTExport('+unittest/xsensBodyData/S01-Trial-Static-1', ...
                                       options);
r = [9.654742,-0.718345,2.264074,-0.047825,0.060515,-0.021490,0.994141,-0.391846,0.300049,0.392628,0.516707,-0.339639,0.680810];
assert(all(dataS.Pelvis.acc(23, :) == r(1:3)));
assert(all(dataS.Pelvis.gyr(23, :) == r(4:6)));
assert(all(dataS.Pelvis.mag(23, :) == r(7:9)));
assert(all(dataS.Pelvis.ori(23, :) == r(10:13)));

%% Test subset data
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataS = mocapdb.XsensBody.loadMTExport('+unittest/xsensBodyData/S01-Trial-Static-1', ...
                                       options);
idx = 1:3:300;
dataSSubset = dataS.getSubset(idx);

for i=1:length(dataSSubset.segList)
    n = dataSSubset.segList{i};
    if sum(size(dataS.(n))) ~= 0
        assert(all(all(table2array(dataSSubset.(n)) == table2array(dataS.(n)(idx, :)))));
    end
end
assert(dataSSubset.nSamples == length(idx));

%% Acceleration normalization
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataS = mocapdb.XsensBody.loadMTExport('+unittest/xsensBodyData/S01-Trial-Static-1', ...
                                       options);
dataSMean = dataS.getMean();
[dataSNorm, dataSBias] = dataS.normalizeAcc();
dataSNorm2 = dataS.removeAccBias(dataSBias);

segList = dataS.getNonemptySeg();
accBiasMean = zeros(1,length(segList));
accBiasNorm = zeros(dataS.nSamples,length(segList));
accBiasNorm2 = zeros(dataS.nSamples,length(segList));
for j = 1:length(segList)
    accBiasMean(:,j) = vecnorm(dataSMean.(segList{j}).acc, 2, 2);
    accBiasNorm(:,j) = vecnorm(dataSNorm.(segList{j}).acc, 2, 2);
    accBiasNorm2(:,j) = vecnorm(dataSNorm2.(segList{j}).acc, 2, 2);
end

tol = 1e-5;
assert(all(all(abs(accBiasNorm-accBiasNorm2)<tol)));
for j = 1:length(segList)
    sname = segList{j};
    assert(all(all(abs(dataS.(sname).acc - ...
        (dataSNorm.(sname).acc + dataSBias.(sname).acc))<tol)));
end