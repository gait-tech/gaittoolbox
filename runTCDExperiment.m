% ======================================================================
%> @brief Run experiment on an instance of the TCD dataset
%> @author Luke Sy
%> 
%> 2018 April 19 changes: Removing support for grlib.est.kf_3_kmus_v2.
%>
%> Setup parameters:
%> - label: data instance name (e.g. s1-acting1 or s2-walking1)
%> - est: filter type to be used.
%>      - ekfv3: grlib.est.kf_3_kmus_v3
%> - accData: acceleration data to be used
%>      - v: vicon
%>      - x: xsens
%> - oriData: orientation data to be used
%>      - v: vicon
%>      - x: xsens
%> - stepDetection: step detection algorithm to be used
%>      - v: fixed variance on tibia vicon data
%>      - x: fixed variance on tibia xsens data
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param fnameV tcd bvh filename or loaded tcdlib.BVHBody 
%>              (e.g. totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh)
%> @param fnameS tcd xsens measurement filename or loaded tcdlib.XsensBody
%>              (e.g. totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors)
%> @param fnameCIB filename of the calibration file from sensor to bone frame
%>                  (e.g. totalcapture/imu/s1/s1_acting1_calib_imu_bone.txt)
%> @param fnameCIR filename of the calibration file from sensor to vicon frame
%>                  (e.g. totalcapture/imu/s1/s1_acting1_calib_imu_ref.txt)
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see
%> details above
%> @param savedir filepath to save .mat output/debug files (optional)
% ======================================================================
function results = runTCDExperiment(fnameV, fnameS, fnameCIB, fnameCIR, ...
                                    name, setups, savedir)
    %% Inputs and Input Check
    validateattributes(fnameV, {'char', 'tcdlib.BVHBody'}, {});
    validateattributes(fnameS, {'char', 'tcdlib.XsensBody'}, {});
    validateattributes(fnameCIB, {'char', 'tcdlib.XsensBody'}, {});
    validateattributes(fnameCIR, {'char', 'tcdlib.XsensBody'}, {});
    if nargin <= 6
        savedir = ''
    end
    
    %% Initialization
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
    
    % Initialize other variables
    fs = 60;
    qV2W = rotm2quat([1 0 0; 0 0 -1; 0 1 0]);

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
    val = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
        'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
    m = containers.Map(key, val);
    for i=1:length(key)
        dataS.(key{i}) = dataS.(key{i})(1:nSamples,:);
        dataV.(val{i}) = dataV.(val{i})(1:nSamples,:)/1000;
    end
    dataV.posUnit = 'm';
    sIdx = 1; eIdx = length(dataV.Hips(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
    
    %% Calculate Orientation
    qOri = {};
    
    % orientation from xsens
    qOri.x = {};
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
    qOri.x.PELV = quatmultiply(qPelvisEst0(sIdx:eIdx, :), qTCD2BM);
    qOri.x.LTIB = quatmultiply(qLankleEst0(sIdx:eIdx, :), qTCD2BM);
    qOri.x.RTIB = quatmultiply(qRankleEst0(sIdx:eIdx, :), qTCD2BM);
    
    % orientation from vicon
    qOri.v.PELV = dataV.qHips(sIdx+1:eIdx+1, :);
    qOri.v.LTIB = dataV.qLeftLeg(sIdx+1:eIdx+1, :);
    qOri.v.RTIB = dataV.qRightLeg(sIdx+1:eIdx+1, :);
    
    %% Position, Velocity, Acceleration
    gfrAcc = {};
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
    gfr_vel_LA_act = diff(dataV.LeftFoot, 1, 1)*fs;
    gfr_vel_LA_act = [gfr_vel_LA_act(1,:); gfr_vel_LA_act];
    gfr_vel_RA_act = diff(dataV.RightFoot, 1, 1)*fs;
    gfr_vel_RA_act = [gfr_vel_RA_act(1,:); gfr_vel_RA_act];
    
        
    x0_pos_MP = MIDPEL_act(sIdx,:);
    x0_pos_LA = dataV.LeftFoot(sIdx,:);
    x0_pos_RA = dataV.RightFoot(sIdx,:);
    x0_vel_MP = gfr_vel_MP_act(sIdx,:);
    x0_vel_LA = gfr_vel_LA_act(sIdx,:);
    x0_vel_RA = gfr_vel_RA_act(sIdx,:);
    
    vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
    
    for i = 1:length(vsigma)
        vLabel = getVLabel('v', vsigma(i));
        gfrAcc.(vLabel) = {};
        gfrAcc.(vLabel).MP = [0 0 0; diff(MIDPEL_act, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).LA = [0 0 0; diff(dataV.LeftFoot, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
        gfrAcc.(vLabel).RA = [0 0 0; diff(dataV.RightFoot, 2, 1)*fs*fs] ...
                             + randn(eIdx,3).*vsigma(i);
    end
    
    % gfrAcc from xsens
    gfrAcc.x = {};
    gfrAcc.x.MP = quatrotate(quatconj(dataS.Pelvis.ori), ...
                            dataS.Pelvis.acc) - [0 0 9.81];
    gfrAcc.x.MP = gfrAcc.x.MP(sIdx:eIdx,:);
%     gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
    gfrAcc.x.LA = quatrotate(quatconj(dataS.L_LowLeg.ori), ...
                            dataS.L_LowLeg.acc) - [0 0 9.81];
    gfrAcc.x.LA = gfrAcc.x.LA(sIdx:eIdx,:);
%     gfr_acc_LA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_LA);
    gfrAcc.x.RA = quatrotate(quatconj(dataS.R_LowLeg.ori), ...
                            dataS.R_LowLeg.acc) - [0 0 9.81];
    gfrAcc.x.RA = gfrAcc.x.RA(sIdx:eIdx,:);
%     gfr_acc_RA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_RA);
    
    % gfrAcc from filtered xsens
    fc = 10;
    [lpf_b, lpf_a] = butter(6, fc/(fs/2));
    gfrAcc.xf.MP = filter(lpf_b, lpf_a, gfrAcc.x.MP);
    gfrAcc.xf.LA = filter(lpf_b, lpf_a, gfrAcc.x.LA);
    gfrAcc.xf.RA = filter(lpf_b, lpf_a, gfrAcc.x.RA);
    
    %% UWB measurements
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
    
    %% Run Experiment
    actBody = dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', ...
                         'xyzColor', {'m', 'y', 'c'}}); 
    actBodyRel = actBody.changeRefFrame('MIDPEL');
            
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
        
        try
            if cs.est == 'ekfv3'
                x0 = [x0_pos_MP x0_vel_MP zeros(1,4) ...
                      x0_pos_LA x0_vel_LA zeros(1,4) ...
                      x0_pos_RA x0_vel_RA zeros(1,4)]';
                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigmaQAccMP', cs.sigmaQAcc, ...
                    'sigmaQAccLA', cs.sigmaQAcc, 'sigmaQAccRA', cs.sigmaQAcc);
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = grlib.est.kf_3_kmus_v3( ...
                    x0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
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
                   gfr_vel_MP_act(sIdx:eIdx,:) qOri.v.PELV...
                   dataV.LeftFoot(sIdx:eIdx,:) ...
                   gfr_vel_LA_act(sIdx:eIdx,:) qOri.v.LTIB...
                   dataV.RightFoot(sIdx:eIdx,:) ...
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
            results0 = estBodyRel.diffRelRMSE(actBodyRel);
        catch
            results0 = actBody.diffRelRMSE(nan);
        end
        
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