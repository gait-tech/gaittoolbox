function results = runNeuRASparse01Experiment(dataS, dataV, ...
                        calibV2W, calibYawFix, ...
                        dataX, revStepDetect, ...
                        name, setups, savedir, startFrame, endFrame, bias)
    % Run experiment on an instance of Vicon and Xsens dataset taken at NeuRA institute
    % 
    % Setup parameters:
    % - label: data instance name (e.g. s1-acting1 or s2-walking1)
    % - est: filter type to be used.
    %      - ckf: pelib.est.kf_3_kmus_v3
    % - accData: acceleration data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - w__x: xsens
    % - oriData: orientation data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - w__x: xsens
    % - stepDetection: step detection algorithm to be used
    %      - false: turn off
    %      - av03: use reviewed step detect file. Dimension similar to dataS
    % - initSrc: source of sensor to body orientation and position init
    %      - w__v: vicon (world frame) (default)
    %      - w__x: xsens
    % - applyMeas: measurement configuration number
    % - applyCstr: constraint configuration number
    % - sigmaQAcc: Q acceleration sigma (variance)
    % - P: initial P matrix
    %
    % :param dataS: loaded mocapdb.XsensBody
    % :param dataV: loaded mocapdb.ViconBody
    % :param calibV2W: quaternion (1 x 4) transforming vicon frame to world frame
    % :param calibYawFix: mocapdb.XsensBody fix ankle sensor yaw orientation offset
    % :param dataX: loaded mocapdb.BVHBody 
    % :param revStepDetect: manually reviewed table if stepL and stepR was detected
    % :param name: name of the experiment
    % :param setups: list of experiment parameters (struct) to be run. see details above
    % :param savedir: filepath to save .mat output/debug files (optional)
    % :param startFrame: frame number at which the algorithm will start
    % :param endFrame: frame number at which the algorithm will end
    % :param bias: pelvis accelerometer bias in sensor frame

    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody', 'numeric'}, {});
    validateattributes(calibV2W, {'numeric'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 9, savedir = ''; end
    if nargin <= 10 || startFrame < 0, startFrame = 100; end
    if nargin <= 11 || endFrame < 0, endFrame = inf; end
    if nargin <= 12
        bias = struct('w__v', zeros(1, 3), 'w__x', zeros(1, 3));
    end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    setupDefault = struct('label', 'ckf', 'est', 'ckf', ...
        'accData', 'v', 'accDataNoise', 0.0, 'oriData', 'v', ...
        'initSrc', 'v', 'stepDetection', 'av01', ...
        'stepDetectWindow', 0.25, 'stepDetectThreshold', 1, ...
        'applyPred', 0, 'applyMeas', 0, 'applyCstr', 0, 'accPelvBias', false, ...
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

%     uwbMeasSigma = 0.1; %m standard deviation
    
    allIdx = {};
    qOri = {};
    gfrAcc = {};
    
    wbodyOri = {};
    bodyAcc = {};
    x0 = {};
    uwbMeas = {};
        
    %% Preprocessing in world frame
    if ~isempty(dataV) && ~isempty(calibV2W)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        W__dataV = dataV.getSubset(1:nSamples).toWorldFrame(calibV2W);
        W__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__v;
        % apply yaw offset to orientation
        W__dataS.L_LowLeg.ori = quatmultiply(calibYawFix.L_LowLeg.ori, W__dataS.L_LowLeg.ori);
        W__dataS.R_LowLeg.ori = quatmultiply(calibYawFix.R_LowLeg.ori, W__dataS.R_LowLeg.ori);
        
        sIdx = max(W__dataV.getStartIndex()+1, startFrame);
        eIdx = min(length(W__dataV.PELV(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__v = idx;
        
        viconCalibSB = W__dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));       
        %% orientation and angular velocity
        % Orientation of body in world frame as obtained from sparse sensor
        qPelvisEst0 = quatmultiply(W__dataS.Pelvis.ori, quatconj(viconCalibSB.Pelvis.ori));
        qLankleEst0 = quatmultiply(W__dataS.L_LowLeg.ori, quatconj(viconCalibSB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(W__dataS.R_LowLeg.ori, quatconj(viconCalibSB.R_LowLeg.ori));
        qOri.w__sv.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.w__sv.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.w__sv.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % https://math.stackexchange.com/questions/2282938/converting-from-quaternion-to-angular-velocity-then-back-to-quaternion
        % Angular velocity of body in body frame as obtained from sparse sensor
        wbodyOri.w__sv.PELV = quatrotate(quatconj(viconCalibSB.Pelvis.ori), W__dataS.Pelvis.gyr);
        wbodyOri.w__sv.PELV = wbodyOri.w__sv.PELV(sIdx:eIdx, :);
        wbodyOri.w__sv.LTIB = quatrotate(quatconj(viconCalibSB.L_LowLeg.ori), W__dataS.L_LowLeg.gyr);
        wbodyOri.w__sv.LTIB = wbodyOri.w__sv.LTIB(sIdx:eIdx, :);
        wbodyOri.w__sv.RTIB = quatrotate(quatconj(viconCalibSB.R_LowLeg.ori), W__dataS.R_LowLeg.gyr);
        wbodyOri.w__sv.RTIB = wbodyOri.w__sv.RTIB(sIdx:eIdx, :);
        
        % Orientation of body in world frame as obtained from vicon input
        qOri.w__v.PELV = W__dataV.qRPV(sIdx+1:eIdx+1, :);
        qOri.w__v.LTIB = W__dataV.qLSK(sIdx+1:eIdx+1, :);
        qOri.w__v.RTIB = W__dataV.qRSK(sIdx+1:eIdx+1, :);
        
        % Angular velocity of body in body frame as obtained from vicon input
        W__viconBody = W__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        angvel = W__viconBody.calcSegAngVel({'qRPV', 'qLSK', 'qRSK'}, 'B');
        wbodyOri.w__v.PELV = angvel.qRPV(sIdx+1:eIdx+1,:);
        wbodyOri.w__v.LTIB = angvel.qLSK(sIdx+1:eIdx+1,:);
        wbodyOri.w__v.RTIB = angvel.qRSK(sIdx+1:eIdx+1,:);
                
        %% position, velocity, acceleration
        vel = W__viconBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
        acc = W__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});

        x0.w__v = [W__viconBody.MIDPEL(sIdx,:) vel.MIDPEL(sIdx,:) zeros(1,4) ...
                   W__viconBody.LTIO(sIdx,:) vel.LTIO(sIdx,:) zeros(1,4) ...
                   W__viconBody.RTIO(sIdx,:) vel.RTIO(sIdx,:) zeros(1,4)]';     
        
        vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
        randnN = size(acc.MIDPEL, 1);
        for i = 1:length(vsigma)
            vLabel = getVLabel('w__v', vsigma(i));
            gfrAcc.(vLabel) = {};
            gfrAcc.(vLabel).MP = acc.MIDPEL + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).MP = gfrAcc.(vLabel).MP(sIdx:eIdx,:);
            gfrAcc.(vLabel).LA = acc.LTIO + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).LA = gfrAcc.(vLabel).LA(sIdx:eIdx,:);
            gfrAcc.(vLabel).RA = acc.RTIO + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).RA = gfrAcc.(vLabel).RA(sIdx:eIdx,:);
        end

        % gfrAcc from sparse
        gfrAcc.w__sv = {};
        gfrAcc.w__sv.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        gfrAcc.w__sv.MP = gfrAcc.w__sv.MP(sIdx:eIdx,:);
        gfrAcc.w__sv.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                                W__dataS.L_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__sv.LA = gfrAcc.w__sv.LA(sIdx:eIdx,:);
        gfrAcc.w__sv.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                                W__dataS.R_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__sv.RA = gfrAcc.w__sv.RA(sIdx:eIdx,:);
        
        % gfrAcc from filtered sparse
        fc = 10;
        [lpf_b, lpf_a] = butter(6, fc/(fs/2));
        gfrAcc.w__sfv.MP = filter(lpf_b, lpf_a, gfrAcc.w__sv.MP);
        gfrAcc.w__sfv.LA = filter(lpf_b, lpf_a, gfrAcc.w__sv.LA);
        gfrAcc.w__sfv.RA = filter(lpf_b, lpf_a, gfrAcc.w__sv.RA);
    
        %% body acceleration
        bodyAcc.w__v.MP = quatrotate(qOri.w__v.PELV, gfrAcc.w__v.MP + [0 0 9.81]);
        bodyAcc.w__v.LA = quatrotate(qOri.w__v.LTIB, gfrAcc.w__v.LA + [0 0 9.81]);
        bodyAcc.w__v.RA = quatrotate(qOri.w__v.RTIB, gfrAcc.w__v.RA + [0 0 9.81]);

        bodyAcc.w__sv.MP = quatrotate(quatconj(viconCalibSB.Pelvis.ori), W__dataS.Pelvis.acc);
        bodyAcc.w__sv.MP = bodyAcc.w__sv.MP(sIdx:eIdx, :);
        bodyAcc.w__sv.LA = quatrotate(quatconj(viconCalibSB.L_LowLeg.ori), W__dataS.L_LowLeg.acc);
        bodyAcc.w__sv.LA = bodyAcc.w__sv.LA(sIdx:eIdx, :);
        bodyAcc.w__sv.RA = quatrotate(quatconj(viconCalibSB.R_LowLeg.ori), W__dataS.R_LowLeg.acc);
        bodyAcc.w__sv.RA = bodyAcc.w__sv.RA(sIdx:eIdx, :);
    
        % debug purposes
        W__viconBody = W__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__viconBody = W__viconBody.changeRefFrame('MIDPEL');
    end
    
    if ~isempty(dataX)
        nSamples = min(dataX.nSamples, dataS.nSamples);
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        
        dataX = dataX.toWorldFrame(qXsensV2W);
        W__dataX = dataX.getSubset(1:nSamples);
        W__dataX.changePosUnit('m', true);
        
        sIdx = startFrame;
        eIdx = min(length(W__dataX.Hips(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__x = idx;
        
        % debug purposes
        W__xsensBody = W__dataX.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
    end    
    
    %% Save processing
    if ~strcmp(savedir, '')
        if ~isempty(dataX)
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', 'W__xsensBody', ...
                 'gfrAcc', 'qOri', 'bodyAcc', 'wbodyOri', ...
                 'x0', 'uwbMeas', 'allIdx')
        else
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', ...
                 'gfrAcc', 'qOri', 'bodyAcc', 'wbodyOri', ...
                 'x0', 'uwbMeas', 'allIdx')
        end
    end
            
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        idx = allIdx.(cs.initSrc);
        sIdx = idx(1);
        eIdx = idx(end);
        idx0 = 1:(eIdx-sIdx+1);
        
        if cs.accData(end) == 'v'
            csGfrAcc = gfrAcc.(getVLabel(cs.accData, cs.accDataNoise));
            csBodyAcc = bodyAcc.(getVLabel(cs.accData, cs.accDataNoise));
        elseif ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
           strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
            csGfrAcc = gfrAcc.(strcat(cs.accData, cs.initSrc(end)));
            csBodyAcc = bodyAcc.(strcat(cs.accData, cs.initSrc(end)));
        else
            csGfrAcc = gfrAcc.(cs.accData);
            csBodyAcc = bodyAcc.(cs.accData);
        end
        
        if strcmp(cs.oriData, 'w__s') || strcmp(cs.oriData, 'v__s')
            csQOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
            csBodyWOri = wbodyOri.(strcat(cs.oriData, cs.initSrc(end)));
        else
            csQOri = qOri.(cs.oriData);
            csBodyWOri = wbodyOri.(cs.oriData);
        end
        
        if ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
           strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
            % init velocity adjustment so pos at t=1 is equal to
            % vicon/xsens body
            csx0 = x0.(cs.initSrc);
            dt = 1.0/fs;
            csx0(4:6, :) = (csx0(4:6, :)*dt - 0.5*csGfrAcc.MP(1,:)'*dt^2)/dt;
            csx0(14:16, :) = (csx0(14:16, :)*dt - 0.5*csGfrAcc.LA(1,:)'*dt^2)/dt;
            csx0(24:26, :) = (csx0(24:26, :)*dt - 0.5*csGfrAcc.RA(1,:)'*dt^2)/dt;
        else
            csx0 = x0.(cs.initSrc);
        end
        
        if strcmp(cs.initSrc, 'w__v')
            csActBody = W__viconBody;
            csActBodyRel = PV__viconBody;
            d_pelvis = norm(W__dataV.RFEP(sIdx,:) - W__dataV.LFEP(sIdx,:));
            d_rfemur = norm(W__dataV.RFEP(sIdx,:) - W__dataV.RFEO(sIdx,:));
            d_lfemur = norm(W__dataV.LFEP(sIdx,:) - W__dataV.LFEO(sIdx,:));
            d_rtibia = norm(W__dataV.RFEO(sIdx,:) - W__dataV.RTIO(sIdx,:));
            d_ltibia = norm(W__dataV.LFEO(sIdx,:) - W__dataV.LTIO(sIdx,:));
        else
            csActBody = W__xsensBody;
            csActBodyRel = PV__xsensBody;
            d_pelvis = norm(W__dataX.RightUpLeg(sIdx,:) - W__dataX.LeftUpLeg(sIdx,:));
            d_rfemur = norm(W__dataX.RightUpLeg(sIdx,:) - W__dataX.RightLeg(sIdx,:));
            d_lfemur = norm(W__dataX.LeftUpLeg(sIdx,:) - W__dataX.LeftLeg(sIdx,:));
            d_rtibia = norm(W__dataX.RightLeg(sIdx,:) - W__dataX.RightFoot(sIdx,:));
            d_ltibia = norm(W__dataX.LeftLeg(sIdx,:) - W__dataX.LeftFoot(sIdx,:));
        end
        
        % step detection
        if strcmp(cs.stepDetection, 'av03')
            csNSamples = size(csGfrAcc.MP, 1);
            bIsStatMP = false(csNSamples, 1);
            bIsStatLA = revStepDetect.stepL(idx);
            bIsStatRA = revStepDetect.stepR(idx);
        else
            csNSamples = size(csGfrAcc.MP, 1);
            bIsStatMP = false(csNSamples, 1);
            bIsStatLA = false(csNSamples, 1);
            bIsStatRA = false(csNSamples, 1);
        end
        
        % if the init position is negative knee angle, allow it
        alphaLKmin = csActBodyRel.calcJointAnglesLKnee(1);
        alphaLKmin = min(alphaLKmin(2), 0);
        alphaRKmin = csActBodyRel.calcJointAnglesRKnee(1);
        alphaRKmin = min(alphaRKmin(2), 0);
        
%         try
            if strcmp(cs.est, 'ckf')

                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigma2QAccMP', cs.sigma2QAcc, ...
                    'sigma2QAccLA', cs.sigma2QAcc, 'sigma2QAccRA', cs.sigma2QAcc, ...
                    'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);
%                 disp(sprintf('%.2f %.2f', rad2deg(alphaLKmin), rad2deg(alphaRKmin)));
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.ckf_3imus( ...
                    csx0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    v3Options);

                idx1EndIdx = find(any(isnan(x_pos_v2), 2), 1);
                if isempty(idx1EndIdx)
                    idx1 = idx0; 
                else 
                    idx1 = idx0(1):idx0(idx1EndIdx-1); 
                    csActBody = csActBody.getSubset(idx1);
                    csActBodyRel = csActBodyRel.getSubset(idx1);
                    disp('NaN found in estimator result. Only evaluating until the last non NaN value');
                end
                
                estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
                   'xyzColor', {'r', 'g', 'b'}, 'fs', fs, ...
                   'MIDPEL', x_pos_v2(idx1, 1:3), ...
                   'LFEP', t_dat_v2.LFEP(idx1, :), ...
                   'LFEO', t_dat_v2.LFEO(idx1, :), ...
                   'LTIO', x_pos_v2(idx1, 11:13), ...
                   'RFEP', t_dat_v2.RFEP(idx1, :), ...
                   'RFEO', t_dat_v2.RFEO(idx1, :), ...
                   'RTIO', x_pos_v2(idx1, 21:23), ...
                   'qRPV', x_pos_v2(idx1, 7:10), ...
                   'qLTH', t_dat_v2.qLTH(idx1, :), ...
                   'qRTH', t_dat_v2.qRTH(idx1, :), ...
                   'qLSK', x_pos_v2(idx1, 17:20), ...
                   'qRSK', x_pos_v2(idx1, 27:30));
               
                estState = x_pos_v2;
                estState2 = t_dat_v2;
%                 actState = [MIDPEL_vicon(sIdx:eIdx,:) ...
%                    gfr_vel_MP_vicon(sIdx:eIdx,:) qOri.v.PELV...
%                    dataV.LTIO(sIdx:eIdx,:) ...
%                    gfr_vel_LA_vicon(sIdx:eIdx,:) qOri.v.LTIB...
%                    dataV.RTIO(sIdx:eIdx,:) ...
%                    gfr_vel_RA_vicon(sIdx:eIdx,:) qOri.v.RTIB];
            elseif strcmp(cs.est, 'pfv1')
                v3Options = struct('fs', fs, 'applyPred', cs.applyPred, ...
                                   'applyMeas', cs.applyMeas );
                lkangle = csActBody.calcJointAnglesLKnee(1);
                rkangle = csActBody.calcJointAnglesRKnee(1);
                csx02 = [csx0(1:6, 1); csActBody.qRPV(1,:)'; ...
                        csActBody.calcJointAnglesLHip(1)'; ...
                        csActBody.calcJointAnglesRHip(1)'; ...
                        lkangle(2); rkangle(2)];
                vel0 = struct('MP', csx0(4:6, 1)', 'LA', csx0(14:16, 1)', ...
                              'RA', csx0(24:26, 1)');
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.pf_3_kmus_v1( ...
                    csx02, false, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwbMeas.(cs.initSrc), vel0, v3Options);
                
                idx1EndIdx = find(any(isnan(x_pos_v2), 2), 1);
                if isempty(idx1EndIdx)
                    idx1 = idx0; 
                else 
                    idx1 = idx0(1):idx0(idx1EndIdx-1); 
                    csActBody = csActBody.getSubset(idx1);
                    csActBodyRel = csActBodyRel.getSubset(idx1);
                    disp('NaN found in estimator result. Only evaluating until the last non NaN value');
                end
                
                zN = zeros(length(idx1), 1);
                estBody = pelib.grBody.generateBodyFromJointAngles(x_pos_v2(idx1, 1:3), ...
                    x_pos_v2(idx1, 7:10), ...
                    x_pos_v2(idx1, 11:13), x_pos_v2(idx1, 14:16), ...
                    [zN x_pos_v2(idx1, 17) zN], [zN x_pos_v2(idx1, 18) zN], ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia);
                estBody.name = 'est'; 
                estBody.posUnit = 'm'; estBody.oriUnit = 'deg';
                estBody.lnSymbol = '--'; estBody.ptSymbol = 'o';
                estBody.frame = 'world'; estBody.xyzColor = {'r', 'g', 'b'};
                estBody.fs = fs;
                estState = x_pos_v2;
                estState2 = t_dat_v2;
            end
            
            runtime = cputime-t0;
            estBody.nSamples = length(idx1);
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'estState', 'estState2', 'runtime', 'cs')
            end

            estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
            csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    %         results(resultsIdx) = estBody.diffRMSE(csActBody);
            results0a = estBody.diffRMSEandMean(csActBody);
            results0 = estBody2.diffRMSEandMean(csActBody2);
%         catch
%             runtime = cputime-t0;
%             results0 = csActBodyRel.diffRMSEandMean(nan);
%         end
        
        results0.dPosW = results0a.dPos;   
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = runtime;
        results(resultsIdx) = results0;
        fprintf("Index %3d/%3d: Running time: %.4f\n", resultsIdx, setupN, cputime-t0);
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