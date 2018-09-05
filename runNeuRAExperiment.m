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
%>      - x: xsens
%>      - s: sparse
%> - stepDetection: step detection algorithm to be used
%>      - 1: fixed variance on tibia accData
%> - initSrc: source of sensor to body orientation and position init
%>      - v: vicon (default)
%>      - x: xsens
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param fnameV loaded mocapdb.ViconBody 
%> @param fnameX loaded mocapdb.BVHBody 
%> @param fnameS loaded mocapdb.XsensBody
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see
%> details above
%> @param savedir filepath to save .mat output/debug files (optional)
% ======================================================================
function results = runNeuRAExperiment(dataV, dataX, dataS, name, setups, savedir)
    %% Inputs and Input Check
    validateattributes(dataV, {'mocapdb.ViconBody'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody'}, {});
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    
    if nargin <= 4
        savedir = '';
    end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    qViconV2W = rotm2quat(eye(3,3));
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'v', 'accDataNoise', 0.0, 'oriData', 'v', ...
        'initSrc', 'v', ...
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
    dataV = dataV.toWorldFrame(qViconV2W);
    if dataX.nSamples > 0
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        dataX = dataX.toWorldFrame(qXsensV2W);
    end
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
        if dataX.nSamples > 0
            dataX.(val2{i}) = dataX.(val2{i})(1:nSamples,:)/1000;
        end
    end
    dataV.posUnit = 'm';
    if dataX.nSamples > 0
        dataX.posUnit = 'm';
    end
    
    sIdx = dataV.getStartIndex()+1; 
    eIdx = length(dataV.PELV(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
    
    viconBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}}); 
    viconBodyRel = viconBody.changeRefFrame('MIDPEL');
    viconCalibWB = dataS.calcCalibSB(viconBody, 1);
    
    if dataX.nSamples > 0
        xsensBody = dataX.togrBody(idx+1, {'name', 'xsens', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                             'xyzColor', {'m', 'y', 'c'}}); 
        xsensBodyRel = xsensBody.changeRefFrame('MIDPEL');
        xsensCalibWB = dataS.calcCalibSB(xsensBody, 1);
    end
    
    
    %% Calculate Orientation
    qOri = {};
    
    % orientation of body from sparse sensor
    qPelvisEst0 = quatmultiply(dataS.Pelvis.ori, quatconj(viconCalibWB.Pelvis.ori));
    qLankleEst0 = quatmultiply(dataS.L_LowLeg.ori, quatconj(viconCalibWB.L_LowLeg.ori));
    qRankleEst0 = quatmultiply(dataS.R_LowLeg.ori, quatconj(viconCalibWB.R_LowLeg.ori));
    qOri.sv.PELV = qPelvisEst0(sIdx:eIdx, :);
    qOri.sv.LTIB = qLankleEst0(sIdx:eIdx, :);
    qOri.sv.RTIB = qRankleEst0(sIdx:eIdx, :);
    
    if dataX.nSamples > 0
        qPelvisEst0 = quatmultiply(dataS.Pelvis.ori, quatconj(xsensCalibWB.Pelvis.ori));
        qLankleEst0 = quatmultiply(dataS.L_LowLeg.ori, quatconj(xsensCalibWB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(dataS.R_LowLeg.ori, quatconj(xsensCalibWB.R_LowLeg.ori));
        qOri.sx.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.sx.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.sx.RTIB = qRankleEst0(sIdx:eIdx, :);
    end
    
    % orientation of body from vicon
    qOri.v.PELV = dataV.qRPV(sIdx+1:eIdx+1, :);
    qOri.v.LTIB = dataV.qLSK(sIdx+1:eIdx+1, :);
    qOri.v.RTIB = dataV.qRSK(sIdx+1:eIdx+1, :);
    
    % orientation of body from xsens
    if dataX.nSamples > 0
        qOri.x.PELV = dataX.qHips(sIdx:eIdx, :);
        qOri.x.LTIB = dataX.qLeftLeg(sIdx:eIdx, :);
        qOri.x.RTIB = dataX.qRightLeg(sIdx:eIdx, :);
    end
    
    %% Position, Velocity, Acceleration
    gfrAcc = {};
    x0 = {};
    
    % gfrAcc from vicon
    MIDPEL_vicon = [mean([dataV.LFEP(:,1) dataV.RFEP(:,1)], 2),...
                  mean([dataV.LFEP(:,2) dataV.RFEP(:,2)], 2),...
                  mean([dataV.LFEP(:,3) dataV.RFEP(:,3)], 2)];
    gfr_vel_MP_vicon = diff(MIDPEL_vicon, 1, 1)*fs;
    gfr_vel_MP_vicon = [gfr_vel_MP_vicon(1,:); gfr_vel_MP_vicon];
    gfr_vel_LA_vicon = diff(dataV.LTIO, 1, 1)*fs;
    gfr_vel_LA_vicon = [gfr_vel_LA_vicon(1,:); gfr_vel_LA_vicon];
    gfr_vel_RA_vicon = diff(dataV.RTIO, 1, 1)*fs;
    gfr_vel_RA_vicon = [gfr_vel_RA_vicon(1,:); gfr_vel_RA_vicon];
    
    x0.v = [MIDPEL_vicon(sIdx,:) gfr_vel_MP_vicon(sIdx,:) zeros(1,4) ...
          dataV.LTIO(sIdx,:) gfr_vel_LA_vicon(sIdx,:) zeros(1,4) ...
          dataV.RTIO(sIdx,:) gfr_vel_RA_vicon(sIdx,:) zeros(1,4)]';     
      
    vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
    for i = 1:length(vsigma)
        vLabel = getVLabel('v', vsigma(i));
        gfrAcc.(vLabel) = {};
        gfrAcc.(vLabel).MP = [0 0 0; diff(MIDPEL_vicon, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).MP = gfrAcc.(vLabel).MP(sIdx:eIdx,:);
        gfrAcc.(vLabel).LA = [0 0 0; diff(dataV.LTIO, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).LA = gfrAcc.(vLabel).LA(sIdx:eIdx,:);
        gfrAcc.(vLabel).RA = [0 0 0; diff(dataV.RTIO, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).RA = gfrAcc.(vLabel).RA(sIdx:eIdx,:);
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
    if dataX.nSamples > 0
        MIDPEL_xsens = [mean([dataX.LeftUpLeg(:,1) dataX.RightUpLeg(:,1)], 2),...
                      mean([dataX.LeftUpLeg(:,2) dataX.RightUpLeg(:,2)], 2),...
                      mean([dataX.LeftUpLeg(:,3) dataX.RightUpLeg(:,3)], 2)];
        gfr_vel_MP_xsens = diff(MIDPEL_xsens, 1, 1)*fs;
        gfr_vel_MP_xsens = [gfr_vel_MP_xsens(1,:); gfr_vel_MP_xsens];
        gfr_vel_LA_xsens = diff(dataX.LeftFoot, 1, 1)*fs;
        gfr_vel_LA_xsens = [gfr_vel_LA_xsens(1,:); gfr_vel_LA_xsens];
        gfr_vel_RA_xsens = diff(dataX.RightFoot, 1, 1)*fs;
        gfr_vel_RA_xsens = [gfr_vel_RA_xsens(1,:); gfr_vel_RA_xsens];
    
        gfrAcc.x = {};
        gfrAcc.x.MP = [0 0 0; diff(MIDPEL_xsens, 2, 1)*fs*fs];
        gfrAcc.x.MP = gfrAcc.x.MP(sIdx:eIdx,:);
        gfrAcc.x.LA = [0 0 0; diff(dataX.LeftFoot, 2, 1)*fs*fs];
        gfrAcc.x.LA = gfrAcc.x.LA(sIdx:eIdx,:);
        gfrAcc.x.RA = [0 0 0; diff(dataX.RightFoot, 2, 1)*fs*fs];
        gfrAcc.x.RA = gfrAcc.x.RA(sIdx:eIdx,:);
        
        x0.x = [MIDPEL_vicon(sIdx,:) gfr_vel_MP_xsens(sIdx,:) zeros(1,4) ...
              dataX.LeftFoot(sIdx,:) gfr_vel_LA_xsens(sIdx,:) zeros(1,4) ...
              dataX.RightFoot(sIdx,:) gfr_vel_RA_xsens(sIdx,:) zeros(1,4)]';      
    end
    
    %% UWB measurements
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    
    uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_vicon-dataV.LTIO), 2, 2) ...
        + normrnd(0, 0.1, [nSamples, 1]);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_vicon-dataV.RTIO), 2, 2) ...
        + normrnd(0, 0.1, [nSamples, 1]);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataV.RTIO-dataV.LTIO), 2, 2) ...
        + normrnd(0, 0.1, [nSamples, 1]);
    
    %% Save processing
    if ~strcmp(savedir, '')
        if dataX.nSamples > 0
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'xsensBody', 'viconBody', 'gfrAcc', 'qOri', 'x0')
        else
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'viconBody', 'gfrAcc', 'qOri', 'x0')
        end
    end
            
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        csGfrAcc = gfrAcc.(getVLabel(cs.accData, cs.accDataNoise));
        if cs.oriData == 's'
            csQOri = qOri.(strcat(cs.oriData, cs.initSrc));
        else
            csQOri = qOri.(cs.oriData);
        end
        csx0 = x0.(cs.initSrc);
        if cs.initSrc == 'v'
            csActBodyRel = viconBodyRel;
            d_pelvis = norm(dataV.RFEP(sIdx,:) - dataV.LFEP(sIdx,:));
            d_rfemur = norm(dataV.RFEP(sIdx,:) - dataV.RFEO(sIdx,:));
            d_lfemur = norm(dataV.LFEP(sIdx,:) - dataV.LFEO(sIdx,:));
            d_rtibia = norm(dataV.RFEO(sIdx,:) - dataV.RTIO(sIdx,:));
            d_ltibia = norm(dataV.LFEO(sIdx,:) - dataV.LTIO(sIdx,:));
        else
            csActBodyRel = xsensBodyRel;
            d_pelvis = norm(dataX.RightUpLeg(sIdx,:) - dataX.LeftUpLeg(sIdx,:));
            d_rfemur = norm(dataX.RightUpLeg(sIdx,:) - dataX.RightLeg(sIdx,:));
            d_lfemur = norm(dataX.LeftUpLeg(sIdx,:) - dataX.LeftLeg(sIdx,:));
            d_rtibia = norm(dataX.RightLeg(sIdx,:) - dataX.RightFoot(sIdx,:));
            d_ltibia = norm(dataX.LeftLeg(sIdx,:) - dataX.LeftFoot(sIdx,:));
        end
        
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
                
                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigmaQAccMP', cs.sigmaQAcc, ...
                    'sigmaQAccLA', cs.sigmaQAcc, 'sigmaQAccRA', cs.sigmaQAcc);
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.kf_3_kmus_v3( ...
                    csx0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
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
%                 actState = [MIDPEL_vicon(sIdx:eIdx,:) ...
%                    gfr_vel_MP_vicon(sIdx:eIdx,:) qOri.v.PELV...
%                    dataV.LTIO(sIdx:eIdx,:) ...
%                    gfr_vel_LA_vicon(sIdx:eIdx,:) qOri.v.LTIB...
%                    dataV.RTIO(sIdx:eIdx,:) ...
%                    gfr_vel_RA_vicon(sIdx:eIdx,:) qOri.v.RTIB];
%                 estAcc = [csGfrAcc.MP(sIdx:eIdx,:) ...
%                           csGfrAcc.LA(sIdx:eIdx,:) ...
%                           csGfrAcc.RA(sIdx:eIdx,:)];
%                 actAcc = [gfrAcc.v.MP(sIdx:eIdx,:) ...
%                           gfrAcc.v.LA(sIdx:eIdx,:) ...
%                           gfrAcc.v.RA(sIdx:eIdx,:)];
            end
            
            runtime = cputime-t0;
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'estState', 'estState2', runtime)
            end
    %         results(resultsIdx) = estBody.diffRMSE(csActBody);
            results0 = estBodyRel.diffRMSE(csActBodyRel);
%         catch
%             results0 = csActBodyRel.diffRMSE(nan);
%         end
        
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = runtime;
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