%% Test grBody togrBody
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
body = dataV.togrBody(1:dataV.nSamples, {});

assert(all(body.LFEP(28, :) == [29.8394000000000,1338.03000000000,804.732000000000]));
assert(all(body.RFEO(14, :) == [-145.434000000000,1380.89000000000,456.709000000000]));
assert(all(body.LTIO(4, :) == [25.3469000000000,1378.58000000000,88.7184000000000]));

%% Test subset data
idx = 1:3:300;
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
dataV = dataV.togrBody(1:dataV.nSamples, {});
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
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
dataV = dataV.togrBody(1:dataV.nSamples, {});
dataV2 = dataV.changePosUnit('m');

for i=1:length(dataV.posList)
    n = dataV.posList{i};
    if sum(size(dataV.(n))) ~= 0
        assert(all(all((dataV2.(n)-dataV.(n)/1000)<1e-10)));
    end
end
assert(dataV2.posUnit == 'm');

%% Test joint velocity and acceleration
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
dataV = dataV.togrBody(1:dataV.nSamples, {});

vel = dataV.calcJointVel();
acc = dataV.calcJointAcc();

fs = dataV.fs;

gfr_vel_MP = diff(dataV.MIDPEL, 1, 1)*fs;
gfr_vel_MP = [gfr_vel_MP(1,:); gfr_vel_MP];
gfr_vel_LA = diff(dataV.LTIO, 1, 1)*fs;
gfr_vel_LA = [gfr_vel_LA(1,:); gfr_vel_LA];
gfr_vel_RA = diff(dataV.RTIO, 1, 1)*fs;
gfr_vel_RA = [gfr_vel_RA(1,:); gfr_vel_RA];

assert(all(all((gfr_vel_MP-vel.MIDPEL) == 0)));
assert(all(all((gfr_vel_LA-vel.LTIO) == 0)));
assert(all(all((gfr_vel_RA-vel.RTIO) == 0)));

gfr_acc_MP = [0 0 0; diff(dataV.MIDPEL, 2, 1)*fs*fs];
gfr_acc_LA = [0 0 0; diff(dataV.LTIO, 2, 1)*fs*fs];
gfr_acc_RA = [0 0 0; diff(dataV.RTIO, 2, 1)*fs*fs];

assert(all(all((gfr_acc_MP-acc.MIDPEL) == 0)));
assert(all(all((gfr_acc_LA-acc.LTIO) == 0)));
assert(all(all((gfr_acc_RA-acc.RTIO) == 0)));