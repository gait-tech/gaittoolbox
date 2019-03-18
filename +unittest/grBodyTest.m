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

%% Test joint angles
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
dataV = dataV.togrBody(1:dataV.nSamples, {});

load('+unittest/grBodyData/S03-Trial-002.mat');
THRESHOLD = 1;
assert(all(all((rad2deg(dataV.calcJointAnglesLHip()) - LHipAngles(:,[2 1 3])) < THRESHOLD)));
assert(all(all((rad2deg(dataV.calcJointAnglesRHip()) - RHipAngles(:,[2 1 3])) < THRESHOLD)));
assert(all(all((rad2deg(dataV.calcJointAnglesLKnee()) - LKneeAngles(:,[2 1 3])) < THRESHOLD)));
assert(all(all((rad2deg(dataV.calcJointAnglesRKnee()) - RKneeAngles(:,[2 1 3])) < THRESHOLD)));

%% Test calcDistR
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
vb = dataV.togrBody(1:dataV.nSamples, {});

LHipAngle = vb.calcJointAnglesLHip(1);
RHipAngle = vb.calcJointAnglesRHip(1);
LKneeAngle = vb.calcJointAnglesLKnee(1);
RKneeAngle = vb.calcJointAnglesRKnee(1);

qLTH = pelib.grBody.calcDistR(vb.qRPV(1,:), LHipAngle.*[-1 -1 -1]);
qRTH = pelib.grBody.calcDistR(vb.qRPV(1,:), RHipAngle.*[1 -1 1]);
qLSK = pelib.grBody.calcDistR(vb.qLTH(1,:), LKneeAngle.*[-1 1 -1]);
qRSK = pelib.grBody.calcDistR(vb.qRTH(1,:), RKneeAngle);

LHipAngle2 = pelib.grBody.calcJointAngles(vb.qRPV(1, :), vb.qLTH(1, :));
LHipAngle2 = LHipAngle2(:, [2 1 3]) .* [-1 -1 -1];
RHipAngle2 = pelib.grBody.calcJointAngles(vb.qRPV(1, :), vb.qRTH(1, :));
RHipAngle2 = RHipAngle2(:, [2 1 3]) .* [1 -1 1];
LKneeAngle2 = pelib.grBody.calcJointAngles(vb.qLTH(1, :), vb.qLSK(1, :));
LKneeAngle2 = LKneeAngle2(:, [2 1 3]) .* [-1 1 -1];
RKneeAngle2 = pelib.grBody.calcJointAngles(vb.qRTH(1, :), vb.qRSK(1, :));
RKneeAngle2 = RKneeAngle2(:, [2 1 3]);
            
THRESHOLD = 1;
assert(all(all(rad2deg(LHipAngle-LHipAngle2) < THRESHOLD)));
assert(all(all(rad2deg(RHipAngle-RHipAngle2) < THRESHOLD)));
assert(all(all(rad2deg(LKneeAngle-LKneeAngle2) < THRESHOLD)));
assert(all(all(rad2deg(RKneeAngle-RKneeAngle2) < THRESHOLD)));

%% Test generateBodyFromJointAngles
dataV = mocapdb.ViconBody.loadViconMat('+unittest/grBodyData/S03-Trial-002.mat');
vb = dataV.togrBody(1:dataV.nSamples, {});

LHipAngle = vb.calcJointAnglesLHip();
RHipAngle = vb.calcJointAnglesRHip();
LKneeAngle = vb.calcJointAnglesLKnee();
RKneeAngle = vb.calcJointAnglesRKnee();
dPelvis = norm(vb.RFEP(1,:) - vb.LFEP(1,:));
dLFemur = norm(vb.LFEP(1,:) - vb.LFEO(1,:));
dRFemur = norm(vb.RFEP(1,:) - vb.RFEO(1,:));
dLTibia = norm(vb.LFEO(1,:) - vb.LTIO(1,:));
dRTibia = norm(vb.RFEO(1,:) - vb.RTIO(1,:));
            
vb2 = pelib.grBody.generateBodyFromJointAngles(vb.MIDPEL, vb.qRPV, ...
    LHipAngle, RHipAngle, LKneeAngle, RKneeAngle, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);


LHipAngle2 = vb.calcJointAnglesLHip();
RHipAngle2 = vb.calcJointAnglesRHip();
LKneeAngle2 = vb.calcJointAnglesLKnee();
RKneeAngle2 = vb.calcJointAnglesRKnee();

THRESHOLD = 1;
assert(all(all(rad2deg(LHipAngle-LHipAngle2) < THRESHOLD)));
assert(all(all(rad2deg(RHipAngle-RHipAngle2) < THRESHOLD)));
assert(all(all(rad2deg(LKneeAngle-LKneeAngle2) < THRESHOLD)));
assert(all(all(rad2deg(RKneeAngle-RKneeAngle2) < THRESHOLD)));