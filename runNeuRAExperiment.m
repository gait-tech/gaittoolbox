% ======================================================================
%> @brief Run experiment on an instance of Vicon and Xsens dataset taken at
%> NeuRA institute
%> @author Luke Sy
%> 
%>
%> Setup parameters:
%> - label: data instance name (e.g. s1-acting1 or s2-walking1)
%> - est: filter type to be used.
%>      - ekfv3: pelib.est.kf_3_kmus_v3
%> - accData: acceleration data to be used
%>      - v: vicon
%>      - s: sparse
%>      - x: xsens
%> - oriData: orientation data to be used
%>      - v: vicon
%>      - s: sparse
%>      - x: xsens
%> - stepDetection: step detection algorithm to be used
%>      - 1: fixed variance on tibia accData
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param fnameV loaded mocapdb.ViconBody 
%>              (e.g. totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh)
%> @param fnameS loaded mocapdb.XsensBody
%>              (e.g. totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors)
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see
%> details above
%> @param savedir filepath to save .mat output/debug files (optional)
% ======================================================================
function results = runNeuRAExperiment(dataV, dataS, dataX, name, setups, savedir)
    %% Inputs and Input Check
    validateattributes(dataV, {'mocapdb.ViconBody'}, {});
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody'}, {});
    
    if nargin <= 4
        savedir = '';
    end
    
    %% Initialization   
    % Initialize other variables
    fs = 100;
    qV2W = rotm2quat(eye(3,3));
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'v', 'accDataNoise', 0.0, 'oriData', 'v', ...
        'stepDetectWindow', 0.25, 'stepDetectThreshold', 1, ...
        'meas', 0, 'cstr', 0, ...
        'sigmaQAcc', 0.5, 'P', 100);
    setupN = length(setups);
    setupDefaultFN = fieldnames(setupDefault);
        
    for i=1:setupN
        for j=1:length(setupDefaultFN)
            if ~isfield(setups{i}, setupDefaultFN{j})
                setups{i}.(setupDefaultFN{j}) = setupDefault.(setupDefaultFN{j});
            end
        end
    end
    
    % convert vicon to world frame
    dataV = dataV.toWorldFrame(qV2W);
    nSamples = min(dataV.nSamples, dataS.nSamples);
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val1 = {'PELV', 'LFEP', 'RFEP',...
        'LFEO', 'RFEO', 'LTIO', 'RTIO'};
    val2 = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
        'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
    for i=1:length(key)
        dataS.(key{i}) = dataS.(key{i})(1:nSamples,:);
        dataV.(val1{i}) = dataV.(val1{i})(1:nSamples,:)/1000;
        dataX.(val2{i}) = dataX.(val2{i})(1:nSamples,:)/1000;
    end
    dataV.posUnit = 'm';
    dataX.posUnit = 'm';
    
    sIdx = 1; eIdx = length(dataV.PELV(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
    
    actBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', ...
                         'xyzColor', {'m', 'y', 'c'}}); 
    actBodyRel = actBody.changeRefFrame('MIDPEL');
    calibIB = dataS.calcCalibSB(actBody);
    
    %% Calculate Orientation
    qOri = {};
    
    % orientation of body from sparse sensor
    qPelvisEst0 = quatmultiply(dataS.Pelvis.ori, quatconj(calibIB.Pelvis.ori));
    qLankleEst0 = quatmultiply(dataS.L_LowLeg.ori, quatconj(calibIB.L_LowLeg.ori));
    qRankleEst0 = quatmultiply(dataS.R_LowLeg.ori, quatconj(calibIB.R_LowLeg.ori));
    qOri.s.PELV = qPelvisEst0(sIdx:eIdx, :);
    qOri.s.LTIB = qLankleEst0(sIdx:eIdx, :);
    qOri.s.RTIB = qRankleEst0(sIdx:eIdx, :);
    
    % orientation of body from vicon
    qOri.v.PELV = dataV.qRPV(sIdx+1:eIdx+1, :);
    qOri.v.LTIB = dataV.qLSK(sIdx+1:eIdx+1, :);
    qOri.v.RTIB = dataV.qRSK(sIdx+1:eIdx+1, :);
    
    % orientation of body from xsens
    qOri.x.PELV = dataX.qHips(sIdx:eIdx, :);
    qOri.x.LTIB = dataX.qLeftLeg(sIdx:eIdx, :);
    qOri.x.RTIB = dataX.qRightLeg(sIdx:eIdx, :);
    
    %% Position, Velocity, Acceleration
    gfrAcc = {};
    
    % gfrAcc from vicon
    MIDPEL_act = [mean([dataV.LFEP(:,1) dataV.RFEP(:,1)], 2),...
                  mean([dataV.LFEP(:,2) dataV.RFEP(:,2)], 2),...
                  mean([dataV.LFEP(:,3) dataV.RFEP(:,3)], 2)];
    gfr_vel_MP_act = diff(MIDPEL_act, 1, 1)*fs;
    gfr_vel_MP_act = [gfr_vel_MP_act(1,:); gfr_vel_MP_act];
    gfr_vel_LA_act = diff(dataV.LTIO, 1, 1)*fs;
    gfr_vel_LA_act = [gfr_vel_LA_act(1,:); gfr_vel_LA_act];
    gfr_vel_RA_act = diff(dataV.RTIO, 1, 1)*fs;
    gfr_vel_RA_act = [gfr_vel_RA_act(1,:); gfr_vel_RA_act];
    
        
    x0_pos_MP = MIDPEL_act(sIdx,:);
    x0_pos_LA = dataV.LTIO(sIdx,:);
    x0_pos_RA = dataV.RTIO(sIdx,:);
    x0_vel_MP = gfr_vel_MP_act(sIdx,:);
    x0_vel_LA = gfr_vel_LA_act(sIdx,:);
    x0_vel_RA = gfr_vel_RA_act(sIdx,:);
    
    vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
    
    for i = 1:length(vsigma)
        vLabel = getVLabel('v', vsigma(i));
        gfrAcc.(vLabel) = {};
        gfrAcc.(vLabel).MP = [0 0 0; diff(MIDPEL_act, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).LA = [0 0 0; diff(dataV.LTIO, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).RA = [0 0 0; diff(dataV.RTIO, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
    end
    
    % gfrAcc from sparse
    gfrAcc.s = {};
    gfrAcc.s.MP = quatrotate(quatconj(dataS.Pelvis.ori), ...
                            dataS.Pelvis.acc) - [0 0 9.81];
    gfrAcc.s.MP = gfrAcc.s.MP(sIdx:eIdx,:);
%     gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
    gfrAcc.s.LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
                            dataS.L_LowLeg.acc) - [0 0 9.81];
    gfrAcc.s.LA = gfrAcc.s.LA(sIdx:eIdx,:);
%     gfr_acc_LA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_LA);
    gfrAcc.s.RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
                            dataS.R_LowLeg.acc) - [0 0 9.81];
    gfrAcc.s.RA = gfrAcc.s.RA(sIdx:eIdx,:);
%     gfr_acc_RA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_RA);
    
    % gfrAcc from filtered sparse
    fc = 10;
    [lpf_b, lpf_a] = butter(6, fc/(fs/2));
    gfrAcc.sf.MP = filter(lpf_b, lpf_a, gfrAcc.s.MP);
    gfrAcc.sf.LA = filter(lpf_b, lpf_a, gfrAcc.s.LA);
    gfrAcc.sf.RA = filter(lpf_b, lpf_a, gfrAcc.s.RA);
    
    % gfrAcc from xsens
    MIDPEL_Xsens = [mean([dataX.LeftUpLeg(:,1) dataX.RightUpLeg(:,1)], 2),...
                  mean([dataX.LeftUpLeg(:,2) dataX.RightUpLeg(:,2)], 2),...
                  mean([dataX.LeftUpLeg(:,3) dataX.RightUpLeg(:,3)], 2)];
    gfrAcc.x = {};
    gfrAcc.x.MP = [0 0 0; diff(MIDPEL_Xsens, 2, 1)*fs*fs];
    gfrAcc.x.MP = gfrAcc.x.MP(sIdx:eIdx,:);
    gfrAcc.x.LA = [0 0 0; diff(dataX.LeftFoot, 2, 1)*fs*fs];
    gfrAcc.x.LA = gfrAcc.x.LA(sIdx:eIdx,:);
    gfrAcc.x.RA = [0 0 0; diff(dataX.RightFoot, 2, 1)*fs*fs];
    gfrAcc.x.RA = gfrAcc.x.RA(sIdx:eIdx,:);
    
    %% UWB measurements
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    
    uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_act-dataV.LTIO), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_act-dataV.RTIO), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.RTIO-dataV.LTIO), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);

    d_pelvis = norm(dataV.RFEP(sIdx,:) - dataV.LFEP(sIdx,:));
    d_rfemur = norm(dataV.RFEP(sIdx,:) - dataV.RFEO(sIdx,:));
    d_lfemur = norm(dataV.LFEP(sIdx,:) - dataV.LFEO(sIdx,:));
    d_rtibia = norm(dataV.RFEO(sIdx,:) - dataV.RTIO(sIdx,:));
    d_ltibia = norm(dataV.LFEO(sIdx,:) - dataV.LTIO(sIdx,:));
    
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        csGfrAcc = gfrAcc.(getVLabel(cs.accData, cs.accDataNoise));
        csQOri = qOri.(cs.oriData);
        
        % step detection
        VAR_WIN  = floor(fs*cs.stepDetectWindow); % NUM_SAMPLES
        ACC_VAR_THRESH = cs.stepDetectThreshold;

        movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
        bIsStatMP = movVarAcc_pelvis < 0;
        movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
        bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
        movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
        bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
        
%         try
            if cs.est == 'ekfv3'
                x0 = [x0_pos_MP x0_vel_MP zeros(1,4) ...
                      x0_pos_LA x0_vel_LA zeros(1,4) ...
                      x0_pos_RA x0_vel_RA zeros(1,4)]';
                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigmaQAccMP', cs.sigmaQAcc, ...
                    'sigmaQAccLA', cs.sigmaQAcc, 'sigmaQAccRA', cs.sigmaQAcc);
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.kf_3_kmus_v3( ...
                    x0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwb_mea, v3Options);
                
                estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
                   'xyzColor', {'r', 'g', 'b'}, 'fs', fs, ...
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
                   gfr_vel_MP_act(sIdx:eIdx,:) qOri.v.PELV...
                   dataV.LTIO(sIdx:eIdx,:) ...
                   gfr_vel_LA_act(sIdx:eIdx,:) qOri.v.LTIB...
                   dataV.RTIO(sIdx:eIdx,:) ...
                   gfr_vel_RA_act(sIdx:eIdx,:) qOri.v.RTIB];
                estAcc = [csGfrAcc.MP(sIdx:eIdx,:) ...
                          csGfrAcc.LA(sIdx:eIdx,:) ...
                          csGfrAcc.RA(sIdx:eIdx,:)];
                actAcc = [gfrAcc.v.MP(sIdx:eIdx,:) ...
                          gfrAcc.v.LA(sIdx:eIdx,:) ...
                          gfrAcc.v.RA(sIdx:eIdx,:)];
            end
            
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'actBody', 'estState', 'actState', ...
                     'estState2', 'estAcc', 'actAcc')
            end
    %         results(resultsIdx) = estBody.diffRMSE(actBody);
            results0 = estBodyRel.diffRMSE(actBodyRel);
%         catch
%             results0 = actBody.diffRMSE(nan);
%         end
        
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = cputime-t0;
        results(resultsIdx) = results0;
        display(sprintf("Index %3d/%3d: Running time: %.4f", resultsIdx, setupN, cputime-t0));
        resultsIdx = resultsIdx + 1;
    end
end

function label = getVLabel(data_source, sigma)
    if data_source == 'v'
        if sigma == 0
            label = 'v';
        else
           label = strrep(sprintf('v%.1f', sigma), '.', '');
        end
    else
        label = data_source;
    end
end