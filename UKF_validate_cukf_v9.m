%% Test cukf_v6 on Simulated Data
%% Import Data
clear all;
close all;
clc;

dataSim = 1;
dataTCD = 2;

data = dataSim;

if data == dataSim
    load('s1-acting1-data.mat')
        %% UKF Setup 
    nSense = 3;
    nStatepSense = 16;
    nMeaspSense = 10;
    nStates=nStatepSense*nSense;      %number of states per sensor = 9: xyz pos, vel, acc
    nMeas = nMeaspSense*nSense; % xyz acc_mp
    %sf = 250; %smaple frequency in Hz
    q=0.1;    %std of process noise
    r=0.4;    %std of measurement noise
    Q=q^2*eye(nStates); % covariance of process
    R=r^2*eye(nMeas);        % covariance of measurement
    sIdx0 = 1;
    
    x0 = [actBody.MIDPEL(sIdx0,:)'; gfr_vel_MP_act(sIdx0,:)'; gfr_acc_MP_filt(sIdx0,:)'; qPelvisAct(sIdx0,:)'; dataS.Pelvis.gyr(sIdx0,:)';...
        actBody.LTIO(sIdx0,:)'; gfr_vel_LA_act(sIdx0,:)'; gfr_acc_LA_filt(sIdx0,:)'; qLankleAct(sIdx0,:)'; dataS.L_LowLeg.gyr(sIdx0,:)';...
        actBody.RTIO(sIdx0,:)'; gfr_vel_RA_act(sIdx0,:)'; gfr_acc_RA_filt(sIdx0,:)'; qRankleAct(sIdx0,:)'; dataS.R_LowLeg.gyr(sIdx0,:)'];
    
    P = eye(length(x0));                         % initial state covraiance
    [N_MP,~] = size(gfr_acc_MP);              % total dynamic steps
    gfr_acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
    isConstr = true;
    
    tStart = tic;
    N_MP = 500; %number rof timesteps to test on KF
    %TODO: return removed estimated values of knee and femur pos/or
    trng = sIdx0:sIdx0+N_MP-1;
    [x_rec, xa_rec, qFEM] = grlib.est.cukf_v9(x0,P,Q,R,N_MP,nMeas,gfr_acc(trng,:),fs, qPelvisEst(trng,:),...
        qLankleEst(trng,:), qRankleEst(trng,:),dataS.Pelvis.gyr(trng,:), dataS.L_LowLeg.gyr(trng,:), dataS.R_LowLeg.gyr(trng,:),...
        d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia,...
        isConstr);
    
    runTime = toc(tStart)
    %save used inintermediate nonlinproj dev.
%     save('nonlinproj_startvar','x0','P','Q','R','N_MP','nMeas','gfr_acc','fs', 'qPelvisEst',...
%         'qLankleEst', 'qRankleEst', 'd_pelvis', 'd_lfemur', 'd_rfemur', 'd_ltibia', 'd_rtibia');
%     
    x_pos_v2 = [x_rec(1:3,:)' x_rec(1+nStatepSense:3+nStatepSense,:)' x_rec(1+2*nStatepSense:3+2*nStatepSense,:)'];
    %xa_rec_v2 = [xa_rec(1:3,:)' xa_rec(4:6,:)' xa_rec(7:9,:)' xa_rec(10:12,:)'];
    xa_rec_v2 = xa_rec';
    x_q_v2 = [x_rec(10:13,:)' x_rec(10+nStatepSense:13+nStatepSense,:)' x_rec(10+2*nStatepSense:13+2*nStatepSense,:)'];
    
    qFEM = qFEM';
    %idx for pos_v2
    MP_pos_Idx = [1:3];
    LA_pos_Idx = [4:6];
    RA_pos_Idx = [7:9];
    
    %idx for x_q_v2
    qMP_Idx = [1:4];
    qLTIB_Idx = [5:8];
    qRTIB_Idx = [9:12];
    
    %idx for xa_rec_v2 (augmented state vec)
    LFEP_pos_Idx = [1:3];
    LFEO_pos_Idx = [4:6];
    RFEP_pos_Idx = [7:9];
    RFEO_pos_Idx = [10:12];
    
    %isolate acc for ploting later
    x_acc_v2 = [x_rec(7:9,:)' x_rec(7+nStatepSense:9+nStatepSense,:)' x_rec(7+2*nStatepSense:9+2*nStatepSense,:)'];
    
    %grlib.viz.UKF_Rel_LPVA_plot
    %% construct bodies
    %use tbl_markers.
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
elseif data == dataTCD
    load('cukf_v6_testfile_txpexp01.mat')
%         %% Swap y and z coordinates of pos data for TCD 
%         %swap y and z col of gfr pos. data of tcd dataset
%         swapyzM = [1 0 0; 0 0 1; 0 1 0];
%         %swp mid pelvis coordinates
%         actBody.MIDPEL = swapyzM*actBody.MIDPEL'
%         actBody.MIDPEL = actBody.MIDPEL';
%         %swap left tibia coordinates
%         actBody.LTIO = swapyzM*actBody.LTIO'
%         actBody.LTIO = actBody.LTIO';
%         %swap right tibia coordinates
%         actBody.RTIO = swapyzM*actBody.RTIO'
%         actBody.RTIO = actBody.RTIO';
%        
%         
    %% UKF Setup 
        
    nSense = 3;
    nStatepSense = 16;
    nMeaspSense = 10;
    nStates=nStatepSense*nSense;      %number of states per sensor = 9: xyz pos, vel, acc
    nMeas = nMeaspSense*nSense; % xyz acc_mp
    %sf = 250; %smaple frequency in Hz
    q=0.1;    %std of process noise
    r=0.4;    %std of measurement noise
    Q=q^2*eye(nStates); % covariance of process
    R=r^2*eye(nMeas);        % covariance of measurement
    sIdx0 = 1;
    
    x0 = [actBody.MIDPEL(sIdx0,:)'; gfr_vel_MP_act(sIdx0,:)'; gfr_acc_MP_filt(sIdx0,:)'; qPelvisAct(sIdx0,:)'; dataS.Pelvis.gyr(sIdx0,:)';...
        actBody.LTIO(sIdx0,:)'; gfr_vel_LA_act(sIdx0,:)'; gfr_acc_LA_filt(sIdx0,:)'; qLankleAct(sIdx0,:)'; dataS.L_LowLeg.gyr(sIdx0,:)';...
        actBody.RTIO(sIdx0,:)'; gfr_vel_RA_act(sIdx0,:)'; gfr_acc_RA_filt(sIdx0,:)'; qRankleAct(sIdx0,:)'; dataS.R_LowLeg.gyr(sIdx0,:)'];
    
    P = eye(length(x0));                         % initial state covraiance
    [N_MP,~] = size(gfr_acc_MP);              % total dynamic steps
    gfr_acc = [gfr_acc_MP, gfr_acc_LA, gfr_acc_RA];
    isConstr = true;
    
    tStart = tic;
    N_MP = 1000; %number of timesteps to test on KF
    %TODO: return removed estimated values of knee and femur pos/or
    trng = sIdx0:sIdx0+N_MP-1;
    [x_rec, xa_rec, qFEM] = grlib.est.cukf_v9(x0,P,Q,R,N_MP,nMeas,gfr_acc(trng,:),fs, qPelvisEst(trng,:),...
        qLankleEst(trng,:), qRankleEst(trng,:),dataS.Pelvis.gyr(trng,:), dataS.L_LowLeg.gyr(trng,:), dataS.R_LowLeg.gyr(trng,:),...
        d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia,...
        isConstr);
    
    runTime = toc(tStart)
    %save used inintermediate nonlinproj dev.
%     save('nonlinproj_startvar','x0','P','Q','R','N_MP','nMeas','gfr_acc','fs', 'qPelvisEst',...
%         'qLankleEst', 'qRankleEst', 'd_pelvis', 'd_lfemur', 'd_rfemur', 'd_ltibia', 'd_rtibia');
%     
    x_pos_v2 = [x_rec(1:3,:)' x_rec(1+nStatepSense:3+nStatepSense,:)' x_rec(1+2*nStatepSense:3+2*nStatepSense,:)'];
    %xa_rec_v2 = [xa_rec(1:3,:)' xa_rec(4:6,:)' xa_rec(7:9,:)' xa_rec(10:12,:)'];
    xa_rec_v2 = xa_rec';
    x_q_v2 = [x_rec(10:13,:)' x_rec(10+nStatepSense:13+nStatepSense,:)' x_rec(10+2*nStatepSense:13+2*nStatepSense,:)'];
    
    qFEM = qFEM';
    %idx for pos_v2
    MP_pos_Idx = [1:3];
    LA_pos_Idx = [4:6];
    RA_pos_Idx = [7:9];
    
    %idx for x_q_v2
    qMP_Idx = [1:4];
    qLTIB_Idx = [5:8];
    qRTIB_Idx = [9:12];
    
    %idx for xa_rec_v2 (augmented state vec)
    LFEP_pos_Idx = [1:3];
    LFEO_pos_Idx = [4:6];
    RFEP_pos_Idx = [7:9];
    RFEO_pos_Idx = [10:12];
    
    %isolate acc for ploting later
    x_acc_v2 = [x_rec(7:9,:)' x_rec(7+nStatepSense:9+nStatepSense,:)' x_rec(7+2*nStatepSense:9+2*nStatepSense,:)'];
    
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

