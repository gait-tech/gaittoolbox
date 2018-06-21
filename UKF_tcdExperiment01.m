%% Reconstruct motion using kf_3_kmus (KF with constraints)
% Steps:
% - calculate orientation (from Xsens or CAHRS)
% - calculate root position by double integration
% - measure angle and position RMSE
% - reconstruct skeleton to feed into visualization
% Sample Run:
% dataV = loadBVHasAnimatedata('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'mm');
% dataS = loadSensors('Data/s1/Acting1_Xsens_AuxFields.sensors');
% animateTrialSample05(dataV, dataS, 60);
% animateTrialSample05('Data/s1/acting1_BlenderZXY_YmZ.bvh', 'Data/s1/Acting1_Xsens_AuxFields.sensors', 60);

% function results = tcdExperiment01(fnameV, fnameS)
%% Inputs and Input Check
clear all;
close all;
clc;
fnameV = 'totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh';
fnameS = 'totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors';
fnameCIB = 'totalcapture/imu/s1/s1_acting1_calib_imu_bone.txt';
fnameCIR = 'totalcapture/imu/s1/s1_acting1_calib_imu_ref.txt';
% Check function input
validateattributes(fnameV, {'char', 'struct'}, {});
validateattributes(fnameS, {'char', 'struct'}, {});

%% Initialization
% Initialize variables and libraries
fs = 60;

% Load calibration data
calibIB = tcdlib.XsensBody.loadCalib(fnameCIB);
calibIR = tcdlib.XsensBody.loadCalib(fnameCIR);

% Load video orientation and position for each body segment
if ischar(fnameV)
    dataV = tcdlib.BVHBody.loadBVHFile(fnameV, 'mm');
else
    dataV = fnameV;
end

% Load sensor orientation and position for each body segment
if ischar(fnameS)
    dataS = tcdlib.XsensBody.loadSensorFile(fnameS);
else
    dataS = fnameS;
end

n = min(dataV.nSamples, dataS.nSamples);
key = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
    'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
val = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
    'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
m = containers.Map(key, val);
for i=1:length(key)
    dataS.(key{i}) = dataS.(key{i})(1:n,:);
    dataV.(val{i}) = dataV.(val{i})(1:n,:)/1000;
    %         dataV.(val{i}) = quatrotate(calibIR.(key{i}).ori, ...
    %             dataV.(val{i})(1:n,:)/1000);
    %         v2 = strcat('q', val{i});
    %         dataV.(v2) = quatmultiply(quatconj(calibIR.(key{i}).ori), ...
    %             dataV.(v2)(1:n,:)/1000);
end
dataV.posUnit = 'm';

MIDPEL_act = [mean([dataV.LeftUpLeg(:,1) dataV.RightUpLeg(:,1)], 2),...
    mean([dataV.LeftUpLeg(:,2) dataV.RightUpLeg(:,2)], 2),...
    mean([dataV.LeftUpLeg(:,3) dataV.RightUpLeg(:,3)], 2)];
%     gfr_vel_MP_act = [gradient(MIDPEL_act(:,1), 1/fs), ...
%                       gradient(MIDPEL_act(:,2), 1/fs), ...
%                       gradient(MIDPEL_act(:,3), 1/fs)];
%     gfr_acc_MP_act = [gradient(gfr_vel_MP_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_MP_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_MP_act(:,3), 1/fs)];
%     gfr_vel_LA_act = [gradient(dataV.LeftFoot(:,1), 1/fs), ...
%                       gradient(dataV.LeftFoot(:,2), 1/fs), ...
%                       gradient(dataV.LeftFoot(:,3), 1/fs)];
%     gfr_acc_LA_act = [gradient(gfr_vel_LA_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_LA_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_LA_act(:,3), 1/fs)];
%     gfr_vel_RA_act = [gradient(dataV.RightFoot(:,1), 1/fs), ...
%                       gradient(dataV.RightFoot(:,2), 1/fs), ...
%                       gradient(dataV.RightFoot(:,3), 1/fs)];
%     gfr_acc_RA_act = [gradient(gfr_vel_RA_act(:,1), 1/fs), ...
%                       gradient(gfr_vel_RA_act(:,2), 1/fs), ...
%                       gradient(gfr_vel_RA_act(:,3), 1/fs)];
gfr_vel_MP_act = diff(MIDPEL_act, 1, 1)*fs;
gfr_vel_MP_act = [gfr_vel_MP_act(1,:); gfr_vel_MP_act];
gfr_acc_MP_act = [0 0 0; diff(MIDPEL_act, 2, 1)*fs*fs];
gfr_vel_LA_act = diff(dataV.LeftFoot, 1, 1)*fs;
gfr_vel_LA_act = [gfr_vel_LA_act(1,:); gfr_vel_LA_act];
gfr_acc_LA_act = [0 0 0; diff(dataV.LeftFoot, 2, 1)*fs*fs];
gfr_vel_RA_act = diff(dataV.RightFoot, 1, 1)*fs;
gfr_vel_RA_act = [gfr_vel_RA_act(1,:); gfr_vel_RA_act];
gfr_acc_RA_act = [0 0 0; diff(dataV.RightFoot, 2, 1)*fs*fs];

sIdx = 1;
x0_pos_MP = MIDPEL_act(sIdx,:);
x0_pos_LA = dataV.LeftFoot(sIdx,:);
x0_pos_RA = dataV.RightFoot(sIdx,:);
x0_vel_MP = gfr_vel_MP_act(sIdx,:);
x0_vel_LA = gfr_vel_LA_act(sIdx,:);
x0_vel_RA = gfr_vel_RA_act(sIdx,:);
sigma_acc = 0.5;

WIN_SECS = 0.25;
VAR_WIN  = floor(fs*WIN_SECS); % NUM_SAMPLES
ACC_VAR_THRESH = 1;

movVarAcc_pelvis = movingvar(sqrt( sum(dataS.Pelvis.acc.^2,2)) ,VAR_WIN);
bIsStatMP = movVarAcc_pelvis < 0;
movVarAcc_lankle = movingvar(sqrt( sum(dataS.L_LowLeg.acc.^2,2)) ,VAR_WIN);
bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
movVarAcc_rankle = movingvar(sqrt( sum(dataS.R_LowLeg.acc.^2,2)) ,VAR_WIN);
bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;

movVarAcc_pelvis_act = movingvar(sqrt( sum(gfr_acc_MP_act.^2,2)) ,VAR_WIN);
bIsStatMP_act = movVarAcc_pelvis_act < 0;
movVarAcc_lankle_act = movingvar(sqrt( sum(gfr_acc_LA_act.^2,2)) ,VAR_WIN);
bIsStatLA_act = movVarAcc_lankle_act < ACC_VAR_THRESH;
movVarAcc_rankle_act = movingvar(sqrt( sum(gfr_acc_RA_act.^2,2)) ,VAR_WIN);
bIsStatRA_act = movVarAcc_rankle_act < ACC_VAR_THRESH;

% orientation of body in Vicon frame
qPelvisEst0 = quatmultiply((calibIR.Pelvis.ori), ...
    quatmultiply(dataS.Pelvis.ori, quatconj(calibIB.Pelvis.ori)));
qLankleEst0 = quatmultiply((calibIR.L_LowLeg.ori), ...
    quatmultiply(dataS.L_LowLeg.ori, quatconj(calibIB.L_LowLeg.ori)));
qRankleEst0 = quatmultiply((calibIR.R_LowLeg.ori), ...
    quatmultiply(dataS.R_LowLeg.ori, quatconj(calibIB.R_LowLeg.ori)));
% orientation of body in World frame
%     qPelvisEst = quatmultiply(dataS.Pelvis.ori, ...
%                               quatconj(calibIB.Pelvis.ori));
%     qLankleEst = quatmultiply(dataS.L_LowLeg.ori, ...
%                               quatconj(calibIB.L_LowLeg.ori));
%     qRankleEst = quatmultiply(dataS.R_LowLeg.ori, ...
%                               quatconj(calibIB.R_LowLeg.ori));

% Align body frame convention to biomechanics convention
% FROM: x - right, y - back, z - down
% TO: x - front, y - left, z - up
qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
qPelvisEst = quatmultiply(qPelvisEst0, qTCD2BM);
qLankleEst = quatmultiply(qLankleEst0, qTCD2BM);
qRankleEst = quatmultiply(qRankleEst0, qTCD2BM);
qPelvisAct = quatmultiply(dataV.qHips, qTCD2BM);
qLankleAct = quatmultiply(dataV.qLeftLeg, qTCD2BM);
qRankleAct = quatmultiply(dataV.qRightLeg, qTCD2BM);

gfr_acc_MP = quatrotate(quatconj(dataS.Pelvis.ori), ...
    dataS.Pelvis.acc) - [0 0 9.81];
gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
gfr_acc_LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
    dataS.L_LowLeg.acc) - [0 0 9.81];
gfr_acc_LA = quatrotate(quatconj(calibIR.L_LowLeg.ori), gfr_acc_LA);
gfr_acc_RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
    dataS.R_LowLeg.acc) - [0 0 9.81];
gfr_acc_RA = quatrotate(quatconj(calibIR.R_LowLeg.ori), gfr_acc_RA);

fc = 10;
[lpf_b, lpf_a] = butter(6, fc/(fs/2));
gfr_acc_MP_filt = filter(lpf_b, lpf_a, gfr_acc_MP);
gfr_acc_LA_filt = filter(lpf_b, lpf_a, gfr_acc_LA);
gfr_acc_RA_filt = filter(lpf_b, lpf_a, gfr_acc_RA);

%% -----------------------------------------------------------------------
%  Simulate uwb measurement by generating pairwise combinations, using the
%  origin of each bone segment as the root point
uwb_mea = struct;
uwb_mea.left_tibia_mid_pelvis = vecnormalize(MIDPEL_act-dataV.LeftFoot);
uwb_mea.mid_pelvis_right_tibia = vecnormalize(MIDPEL_act-dataV.RightFoot);
uwb_mea.left_tibia_right_tibia = vecnormalize(dataV.RightFoot-dataV.LeftFoot);

d_pelvis = norm(dataV.RightUpLeg(sIdx,:) - dataV.LeftUpLeg(sIdx,:));
d_rfemur = norm(dataV.RightUpLeg(sIdx,:) - dataV.RightLeg(sIdx,:));
d_lfemur = norm(dataV.LeftUpLeg(sIdx,:) - dataV.LeftLeg(sIdx,:));
d_rtibia = norm(dataV.RightLeg(sIdx,:) - dataV.RightFoot(sIdx,:));
d_ltibia = norm(dataV.LeftLeg(sIdx,:) - dataV.LeftFoot(sIdx,:));

% validation initialization
eIdx = length(gfr_acc_MP_act(:,1));
% eIdx = 240;
idx = sIdx:eIdx;
idx0 = 1:(eIdx-sIdx+1);

actBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
    'lnSymbol', '-', 'ptSymbol', '*', ...
    'xyzColor', {'m', 'y', 'c'}});


setupDefault = struct('accData', 'v', 'oriData', 'v', ...
    'zupt', 0, 'uwb', 0, 'hjc', 0, 'fdist', 0, ...
    'kneeangle', 0, 'accbias', 0, ...
    'Qacc', sigma_acc*100, 'P', eye(18, 18)*100);

setups = { ... % base scenarios
    struct('accData', 'v', 'oriData', 'v'), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 0, 'hjc', 0), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 1, 'hjc', 0), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 0, 'hjc', 1), ...
    ... % different data sets
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'v', 'oriData', 'x', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'xf', 'oriData', 'v', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'xf', 'oriData', 'x', 'zupt', 1, 'hjc', 1), ...
    ... % different hjc
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 2), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 2), ...
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 3), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 3), ...
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 4), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 4), ...
    ... % different acc_bias
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 1, ...
    'accbias', 1, 'P', eye(27, 27)*100), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 1, ...
    'accbias', 1, 'P', eye(27, 27)*100), ...
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 1, ...
    'accbias', 2, 'P', eye(27, 27)*100), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 1, ...
    'accbias', 2, 'P', eye(27, 27)*100), ...
    };
setups = { ... % base scenarios
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 1), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 0, 'hjc', 9), ...
    struct('accData', 'v', 'oriData', 'v', 'zupt', 1, 'hjc', 9), ...
    struct('accData', 'v', 'oriData', 'x', 'zupt', 1, 'hjc', 9), ...
    struct('accData', 'x', 'oriData', 'v', 'zupt', 1, 'hjc', 9), ...
    struct('accData', 'x', 'oriData', 'x', 'zupt', 1, 'hjc', 9), ...
    struct('accData', 'xf', 'oriData', 'v', 'zupt', 1, 'hjc', 9), ...
    struct('accData', 'xf', 'oriData', 'x', 'zupt', 1, 'hjc', 9), ...
    };

resultsIdx = 1; clear results;
%for sI=1:length(setups)
sI = 1;
    cs = struct(setupDefault);
    csF = fieldnames(setups{sI});
    for i=1:length(csF)
        cs.(csF{i}) = setups{sI}.(csF{i});
    end
    
    if cs.accData == 'x'
        cs_gfr_acc_MP = gfr_acc_MP(sIdx:end-1,:);
        cs_gfr_acc_LA = gfr_acc_LA(sIdx:end-1,:);
        cs_gfr_acc_RA = gfr_acc_RA(sIdx:end-1,:);
    elseif cs.accData == 'xf'
        cs_gfr_acc_MP = gfr_acc_MP_filt(sIdx:end-1,:);
        cs_gfr_acc_LA = gfr_acc_LA_filt(sIdx:end-1,:);
        cs_gfr_acc_RA = gfr_acc_RA_filt(sIdx:end-1,:);
    else
        cs_gfr_acc_MP = gfr_acc_MP_act(sIdx:end,:);
        cs_gfr_acc_LA = gfr_acc_LA_act(sIdx:end,:);
        cs_gfr_acc_RA = gfr_acc_RA_act(sIdx:end,:);
    end
    
    if cs.oriData == 'x'
        cs_qPelvis = qPelvisEst(sIdx:end,:);
        cs_qLankle = qLankleEst(sIdx:end,:);
        cs_qRankle = qRankleEst(sIdx:end,:);
    else
        cs_qPelvis = qPelvisAct(sIdx+1:end,:);
        cs_qLankle = qLankleAct(sIdx+1:end,:);
        cs_qRankle = qRankleAct(sIdx+1:end,:);
    end
    
    %         [ x_pri_v2, x_pos_v2, t_dat_v2 ] = grlib.est.kf_3_kmus_v2(fs, ...
    %             cs.Qacc, cs.Qacc, cs.Qacc, cs.P, ...
    %             x0_pos_MP, x0_vel_MP, cs_gfr_acc_MP, ...
    %             bIsStatMP_act(sIdx:end,:), cs_qPelvis, ...
    %             x0_pos_LA, x0_vel_LA, cs_gfr_acc_LA, ...
    %             bIsStatLA_act(sIdx:end,:), cs_qLankle, ...
    %             x0_pos_RA, x0_vel_RA, cs_gfr_acc_RA, ...
    %             bIsStatRA_act(sIdx:end,:), cs_qRankle, ...
    %             d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
    %             cs.zupt, cs.uwb, cs.hjc, cs.fdist, cs.kneeangle, cs.accbias);
    
    %% UKF modifications
    save('cukf_v6_testfile_txpexp01')
    runUKF = false;
    if runUKF
    nSense = 3;
    nStates=9*nSense;      %number of states per sensor = 9: xyz pos, vel, acc
    nMeas = 3*nSense; % xyz acc_mp
    %sf = 250; %smaple frequency in Hz
    q=0.8;    %std of process noise
    r=0.9;    %std of measurement noise
    Q=q^2*eye(nStates); % covariance of process
    R=r^2*eye(nMeas);        % covariance of measurement
    sIdx = 3;
    
    x0 = [x0_pos_MP'; x0_vel_MP'; gfr_acc_MP(sIdx,:)';...
        x0_pos_LA'; x0_vel_LA'; gfr_acc_LA(sIdx,:)';...
        x0_pos_RA'; x0_vel_RA'; gfr_acc_RA(sIdx,:)'];
    
    P = eye(length(x0));                         % initial state covraiance
    [N_MP,~] = size(gfr_acc_MP);              % total dynamic steps
    gfr_acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
    isConstr = true;
    
    tStart = tic;
    %N_MP = 50;
    %TODO: return removed estimated values of knee and femur pos/or
    [x_rec, xa_rec, qFEM] = grlib.est.cukf_v5(x0,P,Q,R,N_MP,nMeas,gfr_acc,fs, qPelvisEst,...
        qLankleEst, qRankleEst, d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia,...
        isConstr);
    
    %save used inintermediate nonlinproj dev.
%     save('nonlinproj_startvar','x0','P','Q','R','N_MP','nMeas','gfr_acc','fs', 'qPelvisEst',...
%         'qLankleEst', 'qRankleEst', 'd_pelvis', 'd_lfemur', 'd_rfemur', 'd_ltibia', 'd_rtibia');
%     
    x_pos_v2 = [x_rec(1:3,:)' x_rec(10:12,:)' x_rec(19:21,:)'];
    xa_rec_v2 = [xa_rec(1:3,:)' xa_rec(4:6,:)' xa_rec(7:9,:)' xa_rec(10:12,:)'];
    qFEM = qFEM';
    %idx for pos_v2
    MP_pos_Idx = [1:3];
    LA_pos_Idx = [4:6];
    RA_pos_Idx = [7:9];
    
    LFEP_pos_Idx = [1:3];
    LFEO_pos_Idx = [4:6];
    RFEP_pos_Idx = [7:9];
    RFEO_pos_Idx = [10:12];
    
    %isolate acc for ploting later
    x_acc_v2 = [x_rec(7:9,:)' x_rec(16:18,:)' x_rec(25:27,:)'];
    
    %grlib.viz.UKF_Rel_LPVA_plot
    %% other stuff
    estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
        'lnSymbol', '--', 'ptSymbol', 'o', ...
        'xyzColor', {'r', 'g', 'b'}, ...
        'MIDPEL', x_pos_v2(idx0,MP_pos_Idx), ...
        'LFEP', xa_rec_v2(idx0,LFEP_pos_Idx), ...
        'LFEO', xa_rec_v2(idx0,LFEO_pos_Idx), ...
        'LTIO', x_pos_v2(idx0,LA_pos_Idx), ...
        'RFEP', xa_rec_v2(idx0,RFEP_pos_Idx), ...
        'RFEO', xa_rec_v2(idx0,RFEO_pos_Idx), ...
        'RTIO', x_pos_v2(idx0,RA_pos_Idx), ...
        'qLTH', qFEM(idx0,[1:4]), ...
        'qRTH', qFEM(idx0,[5:8]), ...
        'qRPV', cs_qPelvis(idx0,:), ...
        'qRSK', cs_qRankle(idx0,:), ...
        'qLSK', cs_qLankle(idx0,:));
    results(resultsIdx) = estBody.diffRMSE(actBody);
    resultsIdx = resultsIdx + 1;
    
% end
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'LA')
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'RA')
grlib.viz.UKF_plotAccel('LA',gfr_acc_LA_filt, x_acc_v2(:,LA_pos_Idx), N_MP, fs)
grlib.viz.UKF_plotAccel('RA',gfr_acc_RA_filt, x_acc_v2(:,RA_pos_Idx), N_MP, fs)
grlib.struct2csvstr(results, true)
%gelib.viz.plotPositionDiff(estBody, actBody, cell{'MIDPEL'; 'LTIO'; 'RTIO'})

    estBody.frame = 'vicon';
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    xlabel('x'); ylabel('y'); zlabel('z');
    estBodyLimitsRel = [estBodyRel.xlim() estBodyRel.ylim() estBodyRel.zlim()];
    for i=idx
        clf; grid;
        xlim(estBodyLimitsRel(1:2)); 
        ylim(estBodyLimitsRel(3:4)); 
        zlim(estBodyLimitsRel(5:6));  
        view(40, 30);
        grlib.viz.plotLowerBody(estBodyRel, i);
        pause(1/1000);
    end

    end
%% previous KF iterations
%     [ x_pri_v2, x_pos_v2, t_dat_v2 ] = kf_3_kmus_v2(fs, ...
%         sigma_acc, sigma_acc, sigma_acc, P, ...
%         x0_pos_MP, x0_vel_MP, gfr_acc_MP(2:end,:), bIsStatMP_act, qPelvisAct, ...
%         x0_pos_LA, x0_vel_LA, gfr_acc_LA(2:end,:), bIsStatLA_act, qLankleAct, ...
%         x0_pos_RA, x0_vel_RA, gfr_acc_RA(2:end,:), bIsStatRA_act, qRankleAct, ...
%         d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
%         true, false, true, false, false);

%     [ x_pri_v2, x_pos_v2, t_dat_v2 ] = kf_3_kmus_v2(fs, ...
%         0, 0, 0, P, ...
%         x0_pos_MP, x0_vel_MP, gfr_acc_MP_act, bIsStatMP_act, qPelvisAct(2:end,:), ...
%         x0_pos_LA, x0_vel_LA, gfr_acc_LA_act, bIsStatLA_act, qLankleAct(2:end,:), ...
%         x0_pos_RA, x0_vel_RA, gfr_acc_RA_act, bIsStatRA_act, qRankleAct(2:end,:), ...
%         d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
%         true, false, true, false, false);

%     actState = [MIDPEL_act gfr_vel_MP_act ...
%                 dataV.LeftFoot gfr_vel_LA_act ...
%                 dataV.RightFoot gfr_vel_RA_act];
%
%     %% --------------------------------------------------------------------
%     %  Further Validation
% Static Plots
%     updateFigureContents('Position');
%     grlib.viz.plotPosition({estBody, actBody}, {'MIDPEL', 'LTIO', 'RTIO'});
%
%     updateFigureContents('Animation Freeze');
%     grid; view(0, 90); hold on;
%     for i=idx0(1):30:idx0(end)
%         grlib.viz.plotLowerBody(estBody, i);
%         grlib.viz.plotLowerBody(actBody, i);
%     end
% %
% %     updateFigureContents('GFR Acc Diff');
% %     diff_gfr_acc_MP = gfr_acc_MP(1:end-1,:) - gfr_acc_MP_act;
% %     diff_gfr_acc_LA = gfr_acc_LA(1:end-1,:) - gfr_acc_LA_act;
% %     diff_gfr_acc_RA = gfr_acc_RA(1:end-1,:) - gfr_acc_RA_act;
% %     grlib.viz.plotXYZ(diff_gfr_acc_MP, diff_gfr_acc_LA, diff_gfr_acc_RA);
%
%     % Animation
%     updateFigureContents('Animation');
%     xlabel('x'); ylabel('y'); zlabel('z');
%     estBodyLimits = [estBody.xlim() estBody.ylim() estBody.zlim()];
%     for i=idx
%         clf; grid;
%         xlim(estBodyLimits(1:2));
%         ylim(estBodyLimits(3:4));
%         zlim(estBodyLimits(5:6));
%         view(0, 180);
%         grlib.viz.plotLowerBody(estBody, i);
%         pause(1/1000);
%     end

%     updateFigureContents('Animation');
%     xlabel('x'); ylabel('y'); zlabel('z');
%     actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()];
%     for i=idx
%         clf; grid;
%         xlim(actBodyLimits(1:2));
%         ylim(actBodyLimits(3:4));
%         zlim(actBodyLimits(5:6));
%         view(0, 180);
%         grlib.viz.plotLowerBody(actBody, i);
%         pause(1/1000);
%     end
% end