% ======================================================================
%> @brief Run experiment on an instance of Vicon and Xsens dataset taken at
%> NeuRA institute
%> @author Luke Sy
%> 
%> Setup parameters:
%> - label: data instance name (e.g. s1-acting1 or s2-walking1)
%> - est: filter type to be used.
%>      - ekfv3: pelib.est.kf_3_kmus_v3
%> - accData: acceleration data to be used
%>      - w__v: vicon (world frame)
%>      - w__s: sparse (world frame)
%>      - v__v: vicon (vicon frame)
%>      - v__s: sparse (vicon frame)
%>      - w__x: xsens
%> - oriData: orientation data to be used
%>      - w__v: vicon (world frame)
%>      - w__s: sparse (world frame)
%>      - v__v: vicon (vicon frame)
%>      - v__s: sparse (vicon frame)
%>      - w__x: xsens
%> - stepDetection: step detection algorithm to be used
%>      - 1: fixed variance on tibia accData
%> - initSrc: source of sensor to body orientation and position init
%>      - w__v: vicon (world frame) (default)
%>      - v__v: vicon (vicon frame)
%>      - w__x: xsens
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param fnameV loaded mocapdb.ViconBody
%> @param calibV2W quaternion (1 x 4) transforming vicon frame to world frame
%> @param calibW2V mocapdb.XsensBody transforming each sensor's world frame
%>                 to vicon frame
%> @param fnameX loaded mocapdb.BVHBody 
%> @param fnameS loaded mocapdb.XsensBody
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see
%> details above
%> @param savedir filepath to save .mat output/debug files (optional)
% ======================================================================
function results = runNeuRAExperiment(dataS, ...
                        dataV, calibV2W, calibW2V, dataX, ...
                        name, setups, savedir)
    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody', 'numeric'}, {});
    validateattributes(calibV2W, {'numeric'}, {});
    validateattributes(calibW2V, {'mocapdb.XsensBody', 'numeric'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 7
        savedir = '';
    end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
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

    allIdx = {};
    qOri = {};
    gfrAcc = {};
    x0 = {};
    
    %% Generate vicon based inputs in world frame
    if ~isempty(dataV) & ~isempty(calibV2W)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        W__dataV = dataV.getSubset(1:nSamples).toWorldFrame(calibV2W);
        W__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        
        sIdx = max(W__dataV.getStartIndex()+1, 100);
        eIdx = length(W__dataV.PELV(:,1)) - 1;
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__v = idx;
        
        viconCalibWB = W__dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));       
        %% orientation
        % orientation of body from sparse sensor
        qPelvisEst0 = quatmultiply(W__dataS.Pelvis.ori, quatconj(viconCalibWB.Pelvis.ori));
        qLankleEst0 = quatmultiply(W__dataS.L_LowLeg.ori, quatconj(viconCalibWB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(W__dataS.R_LowLeg.ori, quatconj(viconCalibWB.R_LowLeg.ori));
        qOri.w__sv.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.w__sv.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.w__sv.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % orientation of body from vicon
        qOri.w__v.PELV = W__dataV.qRPV(sIdx+1:eIdx+1, :);
        qOri.w__v.LTIB = W__dataV.qLSK(sIdx+1:eIdx+1, :);
        qOri.w__v.RTIB = W__dataV.qRSK(sIdx+1:eIdx+1, :);
        
        
        %% position, velocity, acceleration
        W__viconBody = W__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});  
        vel = W__viconBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
        acc = W__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});

        x0.w__v = [W__viconBody.MIDPEL(sIdx,:) vel.MIDPEL(sIdx,:) zeros(1,4) ...
                   W__viconBody.LTIO(sIdx,:) vel.LTIO(sIdx,:) zeros(1,4) ...
                   W__viconBody.RTIO(sIdx,:) vel.RTIO(sIdx,:) zeros(1,4)]';     

        vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
        for i = 1:length(vsigma)
            vLabel = getVLabel('w__v', vsigma(i));
            gfrAcc.(vLabel) = {};
            gfrAcc.(vLabel).MP = acc.MIDPEL + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).MP = gfrAcc.(vLabel).MP(sIdx:eIdx,:);
            gfrAcc.(vLabel).LA = acc.LTIO + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).LA = gfrAcc.(vLabel).LA(sIdx:eIdx,:);
            gfrAcc.(vLabel).RA = acc.RTIO + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).RA = gfrAcc.(vLabel).RA(sIdx:eIdx,:);
        end

        % gfrAcc from sparse
        gfrAcc.w__s = {};
        gfrAcc.w__s.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        gfrAcc.w__s.MP = gfrAcc.w__s.MP(sIdx:eIdx,:);
        gfrAcc.w__s.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                                W__dataS.L_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__s.LA = gfrAcc.w__s.LA(sIdx:eIdx,:);
        gfrAcc.w__s.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                                W__dataS.R_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__s.RA = gfrAcc.w__s.RA(sIdx:eIdx,:);

        % gfrAcc from filtered sparse
        fc = 10;
        [lpf_b, lpf_a] = butter(6, fc/(fs/2));
        gfrAcc.w__sf.MP = filter(lpf_b, lpf_a, gfrAcc.w__s.MP);
        gfrAcc.w__sf.LA = filter(lpf_b, lpf_a, gfrAcc.w__s.LA);
        gfrAcc.w__sf.RA = filter(lpf_b, lpf_a, gfrAcc.w__s.RA);
        
        % debug purposes
        W__viconBody = W__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__viconBody = W__viconBody.changeRefFrame('MIDPEL');
    end
    
    if ~isempty(dataV) & ~isempty(calibW2V)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        V__dataV = dataV.getSubset(1:nSamples);
        V__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        V__dataS = W__dataS.toViconFrame(calibW2V);
        
        sIdx = max(V__dataV.getStartIndex()+1, 100);
        eIdx = length(V__dataV.PELV(:,1)) - 1;
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.v__v = idx;
        
        viconCalibWB = V__dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));       
        %% orientation
        % orientation of body from sparse sensor
        qPelvisEst0 = quatmultiply(V__dataS.Pelvis.ori, quatconj(viconCalibWB.Pelvis.ori));
        qLankleEst0 = quatmultiply(V__dataS.L_LowLeg.ori, quatconj(viconCalibWB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(V__dataS.R_LowLeg.ori, quatconj(viconCalibWB.R_LowLeg.ori));
        qOri.v__sv.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.v__sv.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.v__sv.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % orientation of body from vicon
        qOri.v__v.PELV = V__dataV.qRPV(sIdx+1:eIdx+1, :);
        qOri.v__v.LTIB = V__dataV.qLSK(sIdx+1:eIdx+1, :);
        qOri.v__v.RTIB = V__dataV.qRSK(sIdx+1:eIdx+1, :);
        
        
        %% position, velocity, acceleration
        V__viconBody = V__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});  
        vel = V__viconBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
        acc = V__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});

        x0.v__v = [V__viconBody.MIDPEL(sIdx,:) vel.MIDPEL(sIdx,:) zeros(1,4) ...
                   V__viconBody.LTIO(sIdx,:) vel.LTIO(sIdx,:) zeros(1,4) ...
                   V__viconBody.RTIO(sIdx,:) vel.RTIO(sIdx,:) zeros(1,4)]';     

        vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
        for i = 1:length(vsigma)
            vLabel = getVLabel('v__v', vsigma(i));
            gfrAcc.(vLabel) = {};
            gfrAcc.(vLabel).MP = acc.MIDPEL + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).MP = gfrAcc.(vLabel).MP(sIdx:eIdx,:);
            gfrAcc.(vLabel).LA = acc.LTIO + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).LA = gfrAcc.(vLabel).LA(sIdx:eIdx,:);
            gfrAcc.(vLabel).RA = acc.RTIO + randn(eIdx,3).*vsigma(i);
            gfrAcc.(vLabel).RA = gfrAcc.(vLabel).RA(sIdx:eIdx,:);
        end

        % gfrAcc from sparse
        gfrAcc.v__s = {};
        gfrAcc.v__s.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        gfrAcc.v__s.MP = gfrAcc.v__s.MP(sIdx:eIdx,:);
        gfrAcc.v__s.MP = quatrotate(quatconj(calibW2V.Pelvis.ori), gfrAcc.v__s.MP);
        gfrAcc.v__s.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                                W__dataS.L_LowLeg.acc) - [0 0 9.81];
        gfrAcc.v__s.LA = gfrAcc.v__s.LA(sIdx:eIdx,:);
        gfrAcc.v__s.LA = quatrotate(quatconj(calibW2V.L_LowLeg.ori), gfrAcc.v__s.LA);
        gfrAcc.v__s.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                                W__dataS.R_LowLeg.acc) - [0 0 9.81];
        gfrAcc.v__s.RA = gfrAcc.v__s.RA(sIdx:eIdx,:);
        gfrAcc.v__s.RA = quatrotate(quatconj(calibW2V.R_LowLeg.ori), gfrAcc.v__s.RA);
        
        % gfrAcc from filtered sparse
        fc = 10;
        [lpf_b, lpf_a] = butter(6, fc/(fs/2));
        gfrAcc.v__sf.MP = filter(lpf_b, lpf_a, gfrAcc.v__s.MP);
        gfrAcc.v__sf.LA = filter(lpf_b, lpf_a, gfrAcc.v__s.LA);
        gfrAcc.v__sf.RA = filter(lpf_b, lpf_a, gfrAcc.v__s.RA);
        
        % debug purposes
        V__viconBody = V__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        
        % gfrAcc from sparse
%         gfrAcc.w__s = {};
%         gfrAcc.w__s.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
%                                 W__dataS.Pelvis.acc) - [0 0 9.81];
%         gfrAcc.w__s.MP = gfrAcc.s.MP(sIdx:eIdx,:);
%     %     gfr_acc_MP = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_MP);
%         gfrAcc.w__s.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
%                                 W__dataS.L_LowLeg.acc) - [0 0 9.81];
%         gfrAcc.w__s.LA = gfrAcc.s.LA(sIdx:eIdx,:);
%     %     gfr_acc_LA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_LA);
%         gfrAcc.w__s.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
%                                 W__dataS.R_LowLeg.acc) - [0 0 9.81];
%         gfrAcc.w__s.RA = gfrAcc.s.RA(sIdx:eIdx,:);
%     %     gfr_acc_RA = quatrotate(quatconj(calibIR.Pelvis.ori), gfr_acc_RA);
    end
    
    if ~isempty(dataX)
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        dataX = dataX.toWorldFrame(qXsensV2W);
        val2 = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
                'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
        for i=1:length(key)
            W__dataS.(key{i}) = W__dataS.(key{i})(1:nSamples,:);
            W__dataV.(val1{i}) = W__dataV.(val1{i})(1:nSamples,:)/1000;
            if dataX.nSamples > 0
                dataX.(val2{i}) = dataX.(val2{i})(1:nSamples,:)/1000;
            end
        end
        dataX.posUnit = 'm';

        xsensBody = dataX.togrBody(idx+1, {'name', 'xsens', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                             'xyzColor', {'m', 'y', 'c'}}); 
        xsensBodyRel = xsensBody.changeRefFrame('MIDPEL');
        xsensCalibWB = W__dataS.calcCalibSB(xsensBody, sIdx(1));
    
        %% Calculate Orientation
        qPelvisEst0 = quatmultiply(W__dataS.Pelvis.ori, quatconj(xsensCalibWB.Pelvis.ori));
        qLankleEst0 = quatmultiply(W__dataS.L_LowLeg.ori, quatconj(xsensCalibWB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(W__dataS.R_LowLeg.ori, quatconj(xsensCalibWB.R_LowLeg.ori));
        qOri.sx.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.sx.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.sx.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % orientation of body from xsens
        qOri.x.PELV = dataX.qHips(sIdx:eIdx, :);
        qOri.x.LTIB = dataX.qLeftLeg(sIdx:eIdx, :);
        qOri.x.RTIB = dataX.qRightLeg(sIdx:eIdx, :);

        %% Position, Velocity, Acceleration
        % gfrAcc from xsens
    
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
%     
%     uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_vicon-W__dataV.LTIO), 2, 2) ...
%         + normrnd(0, 0.1, [nSamples, 1]);
%     uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_vicon-W__dataV.RTIO), 2, 2) ...
%         + normrnd(0, 0.1, [nSamples, 1]);
%     uwb_mea.left_tibia_right_tibia = vecnorm((W__dataV.RTIO-W__dataV.LTIO), 2, 2) ...
%         + normrnd(0, 0.1, [nSamples, 1]);
    
    %% Save processing
    if ~strcmp(savedir, '')
        if ~isempty(dataX)
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'xsensBody', 'viconBody', 'gfrAcc', 'qOri', 'x0', 'idx')
        else
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', 'V__viconBody', 'gfrAcc', 'qOri', 'x0', 'idx')
        end
    end
            
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        csGfrAcc = gfrAcc.(getVLabel(cs.accData, cs.accDataNoise));
        if (cs.oriData == 'w__s') | (cs.oriData == 'v__s')
            csQOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
        else
            csQOri = qOri.(cs.oriData);
        end
        csx0 = x0.(cs.initSrc);
        
        if cs.initSrc == 'w__v'
            csActBodyRel = PV__viconBody;
            d_pelvis = norm(W__dataV.RFEP(sIdx,:) - W__dataV.LFEP(sIdx,:));
            d_rfemur = norm(W__dataV.RFEP(sIdx,:) - W__dataV.RFEO(sIdx,:));
            d_lfemur = norm(W__dataV.LFEP(sIdx,:) - W__dataV.LFEO(sIdx,:));
            d_rtibia = norm(W__dataV.RFEO(sIdx,:) - W__dataV.RTIO(sIdx,:));
            d_ltibia = norm(W__dataV.LFEO(sIdx,:) - W__dataV.LTIO(sIdx,:));
        elseif cs.initSrc == 'v__v'
            csActBodyRel = PV__viconBody;
            d_pelvis = norm(V__dataV.RFEP(sIdx,:) - V__dataV.LFEP(sIdx,:));
            d_rfemur = norm(V__dataV.RFEP(sIdx,:) - V__dataV.RFEO(sIdx,:));
            d_lfemur = norm(V__dataV.LFEP(sIdx,:) - V__dataV.LFEO(sIdx,:));
            d_rtibia = norm(V__dataV.RFEO(sIdx,:) - V__dataV.RTIO(sIdx,:));
            d_ltibia = norm(V__dataV.LFEO(sIdx,:) - V__dataV.LTIO(sIdx,:));
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
        
        try
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
            end
            
            runtime = cputime-t0;
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'estState', 'estState2', 'runtime', 'cs')
            end
    %         results(resultsIdx) = estBody.diffRMSE(csActBody);
            results0 = estBodyRel.diffRMSE(csActBodyRel);
        catch
            runtime = cputime-t0;
            results0 = csActBodyRel.diffRMSE(nan);
        end
        
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = runtime;
        results(resultsIdx) = results0;
        display(sprintf("Index %3d/%3d: Running time: %.4f", resultsIdx, setupN, cputime-t0));
        resultsIdx = resultsIdx + 1;
    end
end

function label = getVLabel(data_source, sigma)
    if data_source(end) == 'v'
        if sigma == 0
            label = data_source;
        else
            label = strrep(sprintf('%s%.1f', data_source, sigma), '.', '');
        end
    else
        label = data_source;
    end
end