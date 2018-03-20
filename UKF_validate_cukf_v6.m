%% Test cukf_v6 on Simulated Data
%% Import Data
clear all;
close all;
clc;

dataSim = 1;
dataTCD = 2;

data = dataTCD;

if data == dataSim
    load('cukf_v6_testfile_50F0706.mat')
        %% UKF Setup 
    nSense = 3;
    nStates=13*nSense;      %number of states per sensor = 9: xyz pos, vel, acc
    nMeas = 7*nSense; % xyz acc_mp
    %sf = 250; %smaple frequency in Hz
    q=0.8;    %std of process noise
    r=0.8;    %std of measurement noise
    Q=q^2*eye(nStates); % covariance of process
    R=r^2*eye(nMeas);        % covariance of measurement
    sIdx = 1;
    
    x0 = [tbl_markers.PELO(1,:)'; gfrPelvis.vel(1,:)'; gfrPelvis.acc(1,:)'; qPelvis(1,:)';...
        tbl_markers.LTIO(1,:)'; gfrLankle.vel(1,:)'; gfrLankle.acc(1,:)'; qLankle(1,:)';...
        tbl_markers.RTIO(1,:)'; gfrRankle.vel(1,:)'; gfrRankle.acc(1,:)'; qRankle(1,:)'];
    
    P = eye(length(x0));                         % initial state covraiance
    [N_MP,~] = size(gfr_acc_MP);              % total dynamic steps
    gfr_acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
    isConstr = true;
    
    tStart = tic;
    %N_MP = 50; %number rof timesteps to test on KF
    %TODO: return removed estimated values of knee and femur pos/or
    [x_rec, xa_rec, qFEM] = grlib.est.cukf_v6(x0,P,Q,R,N_MP,nMeas,gfr_acc,fs, qPelvisEst,...
        qLankleEst, qRankleEst, d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia,...
        isConstr);
    
    runTime = toc(tStart)
    %save used inintermediate nonlinproj dev.
%     save('nonlinproj_startvar','x0','P','Q','R','N_MP','nMeas','gfr_acc','fs', 'qPelvisEst',...
%         'qLankleEst', 'qRankleEst', 'd_pelvis', 'd_lfemur', 'd_rfemur', 'd_ltibia', 'd_rtibia');
%     
    x_pos_v2 = [x_rec(1:3,:)' x_rec(14:16,:)' x_rec(27:29,:)'];
    xa_rec_v2 = [xa_rec(1:3,:)' xa_rec(4:6,:)' xa_rec(7:9,:)' xa_rec(10:12,:)'];
    x_q_v2 = [x_rec(10:13,:)' x_rec(23:26,:)' x_rec(36:39,:)'];
    
    qFEM = qFEM';
    %idx for pos_v2
    MP_pos_Idx = [1:3];
    LA_pos_Idx = [4:6];
    RA_pos_Idx = [7:9];
    
    qMP_Idx = [1:4];
    qLTIB_Idx = [5:8];
    qRTIB_Idx = [9:12];
    
    %idx for augmented state vec
    LFEP_pos_Idx = [1:3];
    LFEO_pos_Idx = [4:6];
    RFEP_pos_Idx = [7:9];
    RFEO_pos_Idx = [10:12];
    
    %isolate acc for ploting later
    x_acc_v2 = [x_rec(7:9,:)' x_rec(20:22,:)' x_rec(33:35,:)'];
    
    %grlib.viz.UKF_Rel_LPVA_plot
    %% construct bodies
    %use tbl_markers.
tsidx0 = 1:N_MP;
    actBody = grlib.grBody('name', 'act', 'posUnit', 'm', 'oriUnit', 'deg', ...
        'lnSymbol', '-', 'ptSymbol', '*', ...
        'xyzColor', {'m', 'y', 'c'}, ...
        'MIDPEL', tbl_markers.PELO(tsidx0,:), ...
        'LFEP', tbl_markers.LFEP(tsidx0,:), ...
        'LFEO', tbl_markers.LFEO(tsidx0,:), ...
        'LTIO', tbl_markers.LTIO(tsidx0,:), ...
        'RFEP', tbl_markers.RFEP(tsidx0,:), ...
        'RFEO', tbl_markers.RFEO(tsidx0,:), ...
        'RTIO', tbl_markers.RTIO(tsidx0,:), ...
        'qLTH', quatSegment.Left_Femur(tsidx0,:), ...
        'qRTH', quatSegment.Right_Femur(tsidx0,:), ...
        'qRPV', quatSegment.Mid_Pelvis(tsidx0,:), ...
        'qLSK', quatSegment.Left_Tibia(tsidx0,:), ...
        'qRSK', quatSegment.Right_Tibia(tsidx0,:));
    
    
    tsidx0_ = 1:N_MP;
    estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
        'lnSymbol', '--', 'ptSymbol', 'o', ...
        'xyzColor', {'r', 'g', 'b'}, ...
        'MIDPEL', x_pos_v2(tsidx0_,MP_pos_Idx), ...
        'LFEP', xa_rec_v2(tsidx0_,LFEP_pos_Idx), ...
        'LFEO', xa_rec_v2(tsidx0_,LFEO_pos_Idx), ...
        'LTIO', x_pos_v2(tsidx0_,LA_pos_Idx), ...
        'RFEP', xa_rec_v2(tsidx0_,RFEP_pos_Idx), ...
        'RFEO', xa_rec_v2(tsidx0_,RFEO_pos_Idx), ...
        'RTIO', x_pos_v2(tsidx0_,RA_pos_Idx), ...
        'qLTH', qFEM(tsidx0_,[1:4]), ...
        'qRTH', qFEM(tsidx0_,[5:8]), ...
        'qRPV', x_q_v2(tsidx0_,qMP_Idx), ...
        'qLSK', x_q_v2(tsidx0_,qLTIB_Idx), ...
        'qRSK', x_q_v2(tsidx0_,qRTIB_Idx));
    %results(resultsIdx) = estBody.diffRMSE(actBody);
    %resultsIdx = resultsIdx + 1;
    
% end
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'LA')
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'RA')
%grlib.viz.UKF_plotAccel('LA',gfr_acc_LA_filt, x_acc_v2(:,LA_pos_Idx), N_MP, fs)
%grlib.viz.UKF_plotAccel('RA',gfr_acc_RA_filt, x_acc_v2(:,RA_pos_Idx), N_MP, fs)
%grlib.struct2csvstr(results, true)
%gelib.viz.plotPositionDiff(estBody, actBody, cell{'MIDPEL'; 'LTIO'; 'RTIO'})
elseif data == dataTCD
    load('cukf_v6_testfile_txpexp01.mat')
        %% UKF Setup 
    nSense = 3;
    nStates=13*nSense;      %number of states per sensor = 9: xyz pos, vel, acc
    nMeas = 7*nSense; % xyz acc_mp
    %sf = 250; %smaple frequency in Hz
    q=0.8;    %std of process noise
    r=0.8;    %std of measurement noise
    Q=q^2*eye(nStates); % covariance of process
    R=r^2*eye(nMeas);        % covariance of measurement
    sIdx = 1;
    
    x0 = [actBody.MIDPEL(1,:)'; gfr_vel_MP_act(1,:)'; gfr_acc_MP_filt(1,:)'; qPelvisAct(1,:)';...
        actBody.LTIO(1,:)'; gfr_vel_LA_act(1,:)'; gfr_acc_LA_filt(1,:)'; qLankleAct(1,:)';...
        actBody.RTIO(1,:)'; gfr_vel_RA_act(1,:)'; gfr_acc_RA_filt(1,:)'; qRankleAct(1,:)'];
    
    P = eye(length(x0));                         % initial state covraiance
    [N_MP,~] = size(gfr_acc_MP);              % total dynamic steps
    gfr_acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
    isConstr = true;
    
    tStart = tic;
    %N_MP = 50; %number rof timesteps to test on KF
    %TODO: return removed estimated values of knee and femur pos/or
    [x_rec, xa_rec, qFEM] = grlib.est.cukf_v6(x0,P,Q,R,N_MP,nMeas,gfr_acc,fs, qPelvisEst,...
        qLankleEst, qRankleEst, d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia,...
        isConstr);
    
    runTime = toc(tStart)
    %save used inintermediate nonlinproj dev.
%     save('nonlinproj_startvar','x0','P','Q','R','N_MP','nMeas','gfr_acc','fs', 'qPelvisEst',...
%         'qLankleEst', 'qRankleEst', 'd_pelvis', 'd_lfemur', 'd_rfemur', 'd_ltibia', 'd_rtibia');
%     
    x_pos_v2 = [x_rec(1:3,:)' x_rec(14:16,:)' x_rec(27:29,:)'];
    xa_rec_v2 = [xa_rec(1:3,:)' xa_rec(4:6,:)' xa_rec(7:9,:)' xa_rec(10:12,:)'];
    x_q_v2 = [x_rec(10:13,:)' x_rec(23:26,:)' x_rec(36:39,:)'];
    
    qFEM = qFEM';
    %idx for pos_v2
    MP_pos_Idx = [1:3];
    LA_pos_Idx = [4:6];
    RA_pos_Idx = [7:9];
    
    qMP_Idx = [1:4];
    qLTIB_Idx = [5:8];
    qRTIB_Idx = [9:12];
    
    %idx for augmented state vec
    LFEP_pos_Idx = [1:3];
    LFEO_pos_Idx = [4:6];
    RFEP_pos_Idx = [7:9];
    RFEO_pos_Idx = [10:12];
    
    %isolate acc for ploting later
    x_acc_v2 = [x_rec(7:9,:)' x_rec(20:22,:)' x_rec(33:35,:)'];
    
    %grlib.viz.UKF_Rel_LPVA_plot
    %% construct Est Body (actBody imported) 
    tsidx0_ = 1:N_MP;
    estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
        'lnSymbol', '--', 'ptSymbol', 'o', ...
        'xyzColor', {'r', 'g', 'b'}, ...
        'MIDPEL', x_pos_v2(tsidx0_,MP_pos_Idx), ...
        'LFEP', xa_rec_v2(tsidx0_,LFEP_pos_Idx), ...
        'LFEO', xa_rec_v2(tsidx0_,LFEO_pos_Idx), ...
        'LTIO', x_pos_v2(tsidx0_,LA_pos_Idx), ...
        'RFEP', xa_rec_v2(tsidx0_,RFEP_pos_Idx), ...
        'RFEO', xa_rec_v2(tsidx0_,RFEO_pos_Idx), ...
        'RTIO', x_pos_v2(tsidx0_,RA_pos_Idx), ...
        'qLTH', qFEM(tsidx0_,[1:4]), ...
        'qRTH', qFEM(tsidx0_,[5:8]), ...
        'qRPV', x_q_v2(tsidx0_,qMP_Idx), ...
        'qLSK', x_q_v2(tsidx0_,qLTIB_Idx), ...
        'qRSK', x_q_v2(tsidx0_,qRTIB_Idx));
    %results(resultsIdx) = estBody.diffRMSE(actBody);
    %resultsIdx = resultsIdx + 1;
    
% end
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'LA')
grlib.viz.UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,'RA')
%grlib.viz.UKF_plotAccel('LA',gfr_acc_LA_filt, x_acc_v2(:,LA_pos_Idx), N_MP, fs)
%grlib.viz.UKF_plotAccel('RA',gfr_acc_RA_filt, x_acc_v2(:,RA_pos_Idx), N_MP, fs)
%grlib.struct2csvstr(results, true)
%gelib.viz.plotPositionDiff(estBody, actBody, cell{'MIDPEL'; 'LTIO'; 'RTIO'})
else
    error('select data for filter')
end

