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

function results = tcdExperiment01(fnameV, fnameS, fnameCIB, fnameCIR, ...
                                   name, setups, savedir)
    %% Inputs and Input Check
%     fnameV = 'totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh';
%     fnameS = 'totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors';
%     fnameCIB = 'totalcapture/imu/s1/s1_acting1_calib_imu_bone.txt';
%     fnameCIR = 'totalcapture/imu/s1/s1_acting1_calib_imu_ref.txt';
    % Check function input
    validateattributes(fnameV, {'char', 'tcdlib.BVHBody'}, {});
    validateattributes(fnameS, {'char', 'tcdlib.XsensBody'}, {});
    validateattributes(fnameCIB, {'char', 'tcdlib.XsensBody'}, {});
    validateattributes(fnameCIR, {'char', 'tcdlib.XsensBody'}, {});
    if nargin <= 6
        savedir = ''
    end
    
    %% Initialization
    % Initialize variables and libraries
    fs = 60;
    
    % Load calibration data
    if ischar(fnameCIB)
        calibIB = tcdlib.XsensBody.loadCalib(fnameCIB);
    else
        calibIB = fnameCIB;
    end
    if ischar(fnameCIR)
        calibIR = tcdlib.XsensBody.loadCalib(fnameCIR);
    else
        calibIR = fnameCIR;
    end
    
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
    
    qV2W = rotm2quat([1 0 0; 0 0 -1; 0 1 0]);
    dataV = dataV.toWorldFrame(qV2W);
    
    nSamples = min(dataV.nSamples, dataS.nSamples);
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
        'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
    m = containers.Map(key, val);
    for i=1:length(key)
        dataS.(key{i}) = dataS.(key{i})(1:nSamples,:);
        dataV.(val{i}) = dataV.(val{i})(1:nSamples,:)/1000;
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
    qPelvisEst0 = quatmultiply(qV2W, qPelvisEst0);
    qLankleEst0 = quatmultiply(qV2W, qLankleEst0);
    qRankleEst0 = quatmultiply(qV2W, qRankleEst0);

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
%     gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
    gfr_acc_LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
                            dataS.L_LowLeg.acc) - [0 0 9.81];
%     gfr_acc_LA = quatrotate(quatconj(calibIR.L_LowLeg.ori), gfr_acc_LA);
%     gfr_acc_LA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_LA);
    gfr_acc_RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
                            dataS.R_LowLeg.acc) - [0 0 9.81];
%     gfr_acc_RA = quatrotate(quatconj(calibIR.R_LowLeg.ori), gfr_acc_RA);
%     gfr_acc_RA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_RA);
    
    fc = 10;
    [lpf_b, lpf_a] = butter(6, fc/(fs/2));
    gfr_acc_MP_filt = filter(lpf_b, lpf_a, gfr_acc_MP);
    gfr_acc_LA_filt = filter(lpf_b, lpf_a, gfr_acc_LA);
    gfr_acc_RA_filt = filter(lpf_b, lpf_a, gfr_acc_RA);
    
    %% -----------------------------------------------------------------------
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    
    uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_act-dataV.LeftFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_act-dataV.RightFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.RightFoot-dataV.LeftFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);

    d_pelvis = norm(dataV.RightUpLeg(sIdx,:) - dataV.LeftUpLeg(sIdx,:));
    d_rfemur = norm(dataV.RightUpLeg(sIdx,:) - dataV.RightLeg(sIdx,:));
    d_lfemur = norm(dataV.LeftUpLeg(sIdx,:) - dataV.LeftLeg(sIdx,:));
    d_rtibia = norm(dataV.RightLeg(sIdx,:) - dataV.RightFoot(sIdx,:));
    d_ltibia = norm(dataV.LeftLeg(sIdx,:) - dataV.LeftFoot(sIdx,:));
    
    % validation initialization
    eIdx = length(gfr_acc_MP_act(:,1));
    % eIdx = 480;
    idx = sIdx:eIdx;
    idx0 = 1:(eIdx-sIdx+1);
    
    actBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', ...
                         'xyzColor', {'m', 'y', 'c'}}); 
    actBodyRel = actBody.changeRefFrame('MIDPEL');
    
    setupDefault = struct('label', 'ekfv2', 'est', 'ekfv2', ...
        'accData', 'v', 'oriData', 'v', 'zuptData', 'v', ...
        'zupt', 0, 'uwb', 0, 'hjc', 0, 'fdist', 0, ...
        'kneeangle', 0, 'accbias', 0, 'cstr', 0, ...
        'Qacc', sigma_acc, 'P', 100);
        
    resultsIdx = 1; clear results;
    setupN = length(setups);
    for sI=1:setupN
        t0 = cputime;
        
        cs = struct(setupDefault);
        csF = fieldnames(setups{sI});
        for i=1:length(csF)
            cs.(csF{i}) = setups{sI}.(csF{i});
        end
        
        if cs.accData == 'x'
            cs_gfr_acc_MP = gfr_acc_MP(sIdx:eIdx,:);
            cs_gfr_acc_LA = gfr_acc_LA(sIdx:eIdx,:);
            cs_gfr_acc_RA = gfr_acc_RA(sIdx:eIdx,:);
        elseif cs.accData == 'xf'
            cs_gfr_acc_MP = gfr_acc_MP_filt(sIdx:eIdx,:);
            cs_gfr_acc_LA = gfr_acc_LA_filt(sIdx:eIdx,:);
            cs_gfr_acc_RA = gfr_acc_RA_filt(sIdx:eIdx,:);
        else
            cs_gfr_acc_MP = gfr_acc_MP_act(sIdx:eIdx,:);
            cs_gfr_acc_LA = gfr_acc_LA_act(sIdx:eIdx,:);
            cs_gfr_acc_RA = gfr_acc_RA_act(sIdx:eIdx,:);
        end
        
        if cs.oriData == 'x'
            cs_qPelvis = qPelvisEst(sIdx:eIdx,:);
            cs_qLankle = qLankleEst(sIdx:eIdx,:);
            cs_qRankle = qRankleEst(sIdx:eIdx,:);
        else
            cs_qPelvis = qPelvisAct(sIdx+1:eIdx+1,:);
            cs_qLankle = qLankleAct(sIdx+1:eIdx+1,:);
            cs_qRankle = qRankleAct(sIdx+1:eIdx+1,:);
        end
        
        if cs.zuptData == 'x'
            cs_bIsStatMP = bIsStatMP(sIdx:eIdx,:);
            cs_bIsStatLA = bIsStatLA(sIdx:eIdx,:);
            cs_bIsStatRA = bIsStatRA(sIdx:eIdx,:);
        else
            cs_bIsStatMP = bIsStatMP_act(sIdx:eIdx,:);
            cs_bIsStatLA = bIsStatLA_act(sIdx:eIdx,:);
            cs_bIsStatRA = bIsStatRA_act(sIdx:eIdx,:);
        end
        
%         try
            if cs.est == 'ekfv2'
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = grlib.est.kf_3_kmus_v2(fs, ...
                    cs.Qacc, cs.Qacc, cs.Qacc, cs.P, ...
                    x0_pos_MP, x0_vel_MP, cs_gfr_acc_MP, ...
                    bIsStatMP_act(sIdx:end,:), cs_qPelvis, ...
                    x0_pos_LA, x0_vel_LA, cs_gfr_acc_LA, ...
                    bIsStatLA_act(sIdx:end,:), cs_qLankle, ...
                    x0_pos_RA, x0_vel_RA, cs_gfr_acc_RA, ...
                    bIsStatRA_act(sIdx:end,:), cs_qRankle, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
                    cs.zupt, cs.uwb, cs.hjc, cs.fdist, cs.kneeangle, cs.accbias);
                
                estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'vicon', ...
                   'xyzColor', {'r', 'g', 'b'}, ...
                   'MIDPEL', x_pos_v2(idx0, 1:3), ...
                   'LFEP', t_dat_v2.LFEP(idx0,:), ...
                   'LFEO', t_dat_v2.LFEO(idx0,:), ...
                   'LTIO', x_pos_v2(idx0, 7:9), ...
                   'RFEP', t_dat_v2.RFEP(idx0,:), ...
                   'RFEO', t_dat_v2.RFEO(idx0,:), ...
                   'RTIO', x_pos_v2(idx0, 13:15), ...
                   'qRPV', cs_qPelvis(idx0,:), ...
                   'qRTH', t_dat_v2.qRTH(idx0,:), ...
                   'qLTH', t_dat_v2.qLTH(idx0,:), ...
                   'qRSK', cs_qRankle(idx0,:), ...
                   'qLSK', cs_qLankle(idx0,:));
                
                estState = x_pos_v2;
                estState2 = t_dat_v2;
                actState = [MIDPEL_act(sIdx:eIdx,:) gfr_vel_MP_act(sIdx:eIdx,:) ...
                            dataV.LeftFoot(sIdx:eIdx,:) gfr_vel_LA_act(sIdx:eIdx,:) ...
                            dataV.RightFoot(sIdx:eIdx,:) gfr_vel_RA_act(sIdx:eIdx,:)];
                estAcc = [gfr_acc_MP(sIdx:eIdx,:) ...
                          gfr_acc_LA(sIdx:eIdx,:) ...
                          gfr_acc_RA(sIdx:eIdx,:)];
                actAcc = [gfr_acc_MP_act(sIdx:eIdx,:) ...
                          gfr_acc_LA_act(sIdx:eIdx,:) ...
                          gfr_acc_RA_act(sIdx:eIdx,:)];
            elseif cs.est == 'ekfv3'
                x0 = [x0_pos_MP x0_vel_MP zeros(1,4) ...
                      x0_pos_LA x0_vel_LA zeros(1,4) ...
                      x0_pos_RA x0_vel_RA zeros(1,4)]';
                v3Options = struct('fs', 60, 'applyZupt', cs.zupt, ...
                    'applyUwb', cs.uwb, 'applyAccBias', cs.accbias, ...
                    'applyConst', cs.cstr, 'sigmaQAccMP', cs.Qacc, ...
                    'sigmaQAccLA', cs.Qacc, 'sigmaQAccRA', cs.Qacc);
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = grlib.est.kf_3_kmus_v3( ...
                    x0, cs.P, ...
                    cs_gfr_acc_MP, bIsStatMP_act(sIdx:end,:), cs_qPelvis, ...
                    cs_gfr_acc_LA, bIsStatLA_act(sIdx:end,:), cs_qLankle, ...
                    cs_gfr_acc_RA, bIsStatRA_act(sIdx:end,:), cs_qRankle, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwb_mea, v3Options);
                
                estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'vicon', ...
                   'xyzColor', {'r', 'g', 'b'}, ...
                   'MIDPEL', x_pos_v2(idx0, 1:3), ...
                   'LFEP', t_dat_v2.LFEP(idx0, :), ...
                   'LFEO', t_dat_v2.LFEO(idx0, :), ...
                   'LTIO', x_pos_v2(idx0, 11:13), ...
                   'RFEP', t_dat_v2.RFEP(idx0, :), ...
                   'RFEO', t_dat_v2.RFEO(idx0, :), ...
                   'RTIO', x_pos_v2(idx0, 21:23), ...
                   'qRPV', x_pos_v2(idx0, 7:10), ...
                   'qLTH', t_dat_v2.qLTH(idx0, :), ...
                   'qRTH', t_dat_v2.qRTH(idx0, :), ...
                   'qLSK', x_pos_v2(idx0, 17:20), ...
                   'qRSK', x_pos_v2(idx0, 27:30));
               
                estState = x_pos_v2;
                estState2 = t_dat_v2;
                actState = [MIDPEL_act(sIdx:eIdx,:) ...
                   gfr_vel_MP_act(sIdx:eIdx,:) qPelvisAct(sIdx+1:eIdx+1,:)...
                   dataV.LeftFoot(sIdx:eIdx,:) ...
                   gfr_vel_LA_act(sIdx:eIdx,:) qLankleAct(sIdx+1:eIdx+1,:)...
                   dataV.RightFoot(sIdx:eIdx,:) ...
                   gfr_vel_RA_act(sIdx:eIdx,:) qRankleAct(sIdx+1:eIdx+1,:)];
                estAcc = [gfr_acc_MP(sIdx:eIdx,:) ...
                          gfr_acc_LA(sIdx:eIdx,:) ...
                          gfr_acc_RA(sIdx:eIdx,:)];
                actAcc = [gfr_acc_MP_act(sIdx:eIdx,:) ...
                          gfr_acc_LA_act(sIdx:eIdx,:) ...
                          gfr_acc_RA_act(sIdx:eIdx,:)];
            end
            
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'actBody', 'estState', 'actState', ...
                     'estState2', 'estAcc', 'actAcc')
            end
    %         results(resultsIdx) = estBody.diffRMSE(actBody);
            results0 = estBodyRel.diffRelRMSE(actBodyRel);
%         catch
%             results0 = actBody.diffRelRMSE(nan);
%         end
        
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = cputime-t0;
        results(resultsIdx) = results0;
        display(sprintf("Index %3d/%3d: Running time: %.4f", resultsIdx, setupN, cputime-t0));
        resultsIdx = resultsIdx + 1;
    end
    
%     grlib.viz.plotLowerBodySegmentLengthError(estBody, d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia)
%     struct2table(results)
%     grlib.viz.plotPosition({estBodyRel, actBodyRel}, {'LTIO', 'RTIO'})
%     
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
%     grlib.viz.plotStateComparison(t_dat_v2, actState(1:end-1), 7)
% 
%     %% --------------------------------------------------------------------
%     %  Further Validation
end