% ======================================================================
%> @brief Run experiment on an instance of Vicon and Xsens dataset taken at
%> NeuRA institute
%> @author Luke Sy
%> 
%> Setup parameters:
%> - label: data instance name (e.g. s1-acting1 or s2-walking1)
%> - est: filter type to be used.
%>      - ekfv3: pelib.est.kf_3_kmus_v3
%>      - pfv1: pelib.est.pf_3_kmus_v1
%>      - mpfv1: pelib.est.mpf_3_kmus_v1
%>      - lieekfv1: pelib.est.lieekf_3_kmus_v1
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
%>      - false: turn off
%>      - av01: fixed acceleration variance on tibia accData (var = 1)
%>      - av02: fixed acceleration variance on vicon tibia accData (var = 1)
%>      - av03: use reviewed step detect file. Dimension similar to dataS
%> - initSrc: source of sensor to body orientation and position init
%>      - w__v: vicon (world frame) (default)
%>      - v__v: vicon (vicon frame)
%>      - w__x: xsens
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param dataS loaded mocapdb.XsensBody
%> @param dataV loaded mocapdb.ViconBody
%> @param calibV2W quaternion (1 x 4) transforming vicon frame to world frame
%> @param calibYawFix mocapdb.XsensBody fix ankle sensor yaw orientation offset
%> @param calibW2V mocapdb.XsensBody transforming each sensor's world frame
%>                 to vicon frame. Includes yaw realignment calibration.
%> @param dataX loaded mocapdb.BVHBody 
%> @param revStepDetect manually reviewed table if stepL and stepR was detected
%> @param uwbMeasSigma standard dev of uwb gaussian noise in m
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see details above
%> @param savedir filepath to save .mat output/debug files (optional)
%> @param startFrame frame number at which the algorithm will start
%> @param endFrame frame number at which the algorithm will end
%> @param bias pelvis accelerometer bias in sensor frame
% ======================================================================
function results = runNeuRASparse01Experiment(dataS, dataV, ...
                        calibV2W, calibYawFix, calibW2V, ...
                        dataX, revStepDetect, uwbMeasSigma, ...
                        name, setups, savedir, startFrame, endFrame, bias)
    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody', 'numeric'}, {});
    validateattributes(calibV2W, {'numeric'}, {});
    validateattributes(calibW2V, {'mocapdb.XsensBody', 'numeric'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 7, uwbMeasSigma = 0.0; end
    if nargin <= 10, savedir = ''; end
    if nargin <= 11 || startFrame < 0, startFrame = 100; end
    if nargin <= 12 || endFrame < 0, endFrame = inf; end
    if nargin <= 13
        bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                      'w__x', zeros(1, 3));
    end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'v', 'accDataNoise', 0.0, 'oriData', 'v', ...
        'initSrc', 'v', 'stepDetection', 'av01', ...
        'stepDetectWindow', 0.25, 'stepDetectThreshold', 1, ...
        'meas', 0, 'cstr', 0, 'accPelvBias', false, ...
        'sigmaUwbLeg', 1e0, ...
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
    wOri = {};
    gfrAcc = {};
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
        % orientation of body from sparse sensor
        qPelvisEst0 = quatmultiply(W__dataS.Pelvis.ori, quatconj(viconCalibSB.Pelvis.ori));
        qLankleEst0 = quatmultiply(W__dataS.L_LowLeg.ori, quatconj(viconCalibSB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(W__dataS.R_LowLeg.ori, quatconj(viconCalibSB.R_LowLeg.ori));
        qOri.w__sv.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.w__sv.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.w__sv.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % https://math.stackexchange.com/questions/2282938/converting-from-quaternion-to-angular-velocity-then-back-to-quaternion
        wOri.w__sv.PELV = quatrotate(quatconj(viconCalibSB.Pelvis.ori), W__dataS.Pelvis.gyr);
        wOri.w__sv.PELV = wOri.w__sv.PELV(sIdx:eIdx, :);
        wOri.w__sv.LTIB = quatrotate(quatconj(viconCalibSB.L_LowLeg.ori), W__dataS.L_LowLeg.gyr);
        wOri.w__sv.LTIB = wOri.w__sv.LTIB(sIdx:eIdx, :);
        wOri.w__sv.RTIB = quatrotate(quatconj(viconCalibSB.R_LowLeg.ori), W__dataS.R_LowLeg.gyr);
        wOri.w__sv.RTIB = wOri.w__sv.RTIB(sIdx:eIdx, :);
        
        % orientation of body from vicosIdxn
        qOri.w__v.PELV = W__dataV.qRPV(sIdx+1:eIdx+1, :);
        qOri.w__v.LTIB = W__dataV.qLSK(sIdx+1:eIdx+1, :);
        qOri.w__v.RTIB = W__dataV.qRSK(sIdx+1:eIdx+1, :);
        
        wPelvis0 = quatmultiply(quatconj(W__dataV.qRPV(sIdx+1:eIdx+1, :)), ...
                                W__dataV.qRPV([sIdx+2:eIdx+1 eIdx+1], :));
        wLankle0 = quatmultiply(quatconj(W__dataV.qLSK(sIdx+1:eIdx+1, :)), ...
                                W__dataV.qLSK([sIdx+2:eIdx+1 eIdx+1], :));
        wRankle0 = quatmultiply(quatconj(W__dataV.qRSK(sIdx+1:eIdx+1, :)), ...
                                W__dataV.qRSK([sIdx+2:eIdx+1 eIdx+1], :));
        wOri.w__v.PELV = 2*wPelvis0(:,2:4)*fs;
        wOri.w__v.LTIB = 2*wLankle0(:,2:4)*fs;
        wOri.w__v.RTIB = 2*wRankle0(:,2:4)*fs;
        
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
    
        % UWB measurements
        %  Simulate uwb measurement by generating pairwise combinations, using the
        %  origin of each bone segment as the root point
        umBuf = struct;
        umBuf.left_tibia_mid_pelvis = vecnorm((W__viconBody.MIDPEL-W__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_mid_pelvis = umBuf.left_tibia_mid_pelvis(sIdx:eIdx,:);
        umBuf.mid_pelvis_right_tibia = vecnorm((W__viconBody.MIDPEL-W__viconBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.mid_pelvis_right_tibia = umBuf.mid_pelvis_right_tibia(sIdx:eIdx,:);
        umBuf.left_tibia_right_tibia = vecnorm((W__viconBody.RTIO-W__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_right_tibia = umBuf.left_tibia_right_tibia(sIdx:eIdx,:);
        umBuf.LLeg = vecnorm((W__viconBody.LFEP-W__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.LLeg = umBuf.LLeg(sIdx:eIdx,:);
        umBuf.RLeg = vecnorm((W__viconBody.RFEP-W__viconBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.RLeg = umBuf.RLeg(sIdx:eIdx,:);
        uwbMeas.w__v = umBuf;
    
        % debug purposes
        W__viconBody = W__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__viconBody = W__viconBody.changeRefFrame('MIDPEL');
    end
    
    %% Preprocessing in vicon frame
    if ~isempty(dataV) && ~isempty(calibW2V)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        V__dataV = dataV.getSubset(1:nSamples);
        V__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__v;
        V__dataS = W__dataS.toViconFrame(calibW2V);
        
        sIdx = max(V__dataV.getStartIndex()+1, startFrame);
        eIdx = min(length(V__dataV.PELV(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.v__v = idx;
        
        viconCalibSB = V__dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));       
        %% orientation
        % orientation of body from sparse sensor
        qPelvisEst0 = quatmultiply(V__dataS.Pelvis.ori, quatconj(viconCalibSB.Pelvis.ori));
        qLankleEst0 = quatmultiply(V__dataS.L_LowLeg.ori, quatconj(viconCalibSB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(V__dataS.R_LowLeg.ori, quatconj(viconCalibSB.R_LowLeg.ori));
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
        randnN = size(acc.MIDPEL, 1);
        for i = 1:length(vsigma)
            vLabel = getVLabel('v__v', vsigma(i));
            gfrAcc.(vLabel) = {};
            gfrAcc.(vLabel).MP = acc.MIDPEL + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).MP = gfrAcc.(vLabel).MP(sIdx:eIdx,:);
            gfrAcc.(vLabel).LA = acc.LTIO + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).LA = gfrAcc.(vLabel).LA(sIdx:eIdx,:);
            gfrAcc.(vLabel).RA = acc.RTIO + randn(randnN,3).*vsigma(i);
            gfrAcc.(vLabel).RA = gfrAcc.(vLabel).RA(sIdx:eIdx,:);
        end

        % gfrAcc from sparse
        gfrAcc.v__sv = {};
        gfrAcc.v__sv.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        gfrAcc.v__sv.MP = gfrAcc.v__sv.MP(sIdx:eIdx,:);
        gfrAcc.v__sv.MP = quatrotate(quatconj(calibW2V.Pelvis.ori), gfrAcc.v__sv.MP);
        
        gfrAcc.v__sv.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                                W__dataS.L_LowLeg.acc) - [0 0 9.81];
        gfrAcc.v__sv.LA = gfrAcc.v__sv.LA(sIdx:eIdx,:);
        gfrAcc.v__sv.LA = quatrotate(quatconj(calibW2V.L_LowLeg.ori), gfrAcc.v__sv.LA);
        gfrAcc.v__sv.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                                W__dataS.R_LowLeg.acc) - [0 0 9.81];
        gfrAcc.v__sv.RA = gfrAcc.v__sv.RA(sIdx:eIdx,:);
        gfrAcc.v__sv.RA = quatrotate(quatconj(calibW2V.R_LowLeg.ori), gfrAcc.v__sv.RA);
               
        % gfrAcc from filtered sparse
        fc = 10;
        [lpf_b, lpf_a] = butter(6, fc/(fs/2));
        gfrAcc.v__sfv.MP = filter(lpf_b, lpf_a, gfrAcc.v__sv.MP);
        gfrAcc.v__sfv.LA = filter(lpf_b, lpf_a, gfrAcc.v__sv.LA);
        gfrAcc.v__sfv.RA = filter(lpf_b, lpf_a, gfrAcc.v__sv.RA);
        
        % UWB measurements
        %  Simulate uwb measurement by generating pairwise combinations, using the
        %  origin of each bone segment as the root point
        umBuf = struct;
        umBuf.left_tibia_mid_pelvis = vecnorm((V__viconBody.MIDPEL-V__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_mid_pelvis = umBuf.left_tibia_mid_pelvis(sIdx:eIdx,:);
        umBuf.mid_pelvis_right_tibia = vecnorm((V__viconBody.MIDPEL-V__viconBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.mid_pelvis_right_tibia = umBuf.mid_pelvis_right_tibia(sIdx:eIdx,:);
        umBuf.left_tibia_right_tibia = vecnorm((V__viconBody.RTIO-V__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_right_tibia = umBuf.left_tibia_right_tibia(sIdx:eIdx,:);
        umBuf.LLeg = vecnorm((V__viconBody.LFEP-V__viconBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.LLeg = umBuf.LLeg(sIdx:eIdx,:);
        umBuf.RLeg = vecnorm((V__viconBody.RFEP-V__viconBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.RLeg = umBuf.RLeg(sIdx:eIdx,:);
        uwbMeas.v__v = umBuf;
        
        % debug purposes
        V__viconBody = V__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
    end
    
    if ~isempty(dataX)
        nSamples = min(dataX.nSamples, dataS.nSamples);
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        
        dataX = dataX.toWorldFrame(qXsensV2W);
        W__dataX = dataX.getSubset(1:nSamples);
        W__dataX.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples); % .toViconFrame(calibW2V);
        % order is not important as calibW2V fixes only the ankle yaw offset
        W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__x; 
        % apply yaw offset to orientation
        W__dataS.L_LowLeg.ori = quatmultiply(calibYawFix.L_LowLeg.ori, W__dataS.L_LowLeg.ori);
        W__dataS.R_LowLeg.ori = quatmultiply(calibYawFix.R_LowLeg.ori, W__dataS.R_LowLeg.ori);
        
        sIdx = startFrame;
        eIdx = min(length(W__dataX.Hips(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__x = idx;
        xsensCalibSB = W__dataS.calcCalibSB(W__dataX.togrBody(sIdx+1:sIdx+1, {}), sIdx(1)); 
                
        % orientation of body from sparse sensor
        qPelvisEst0 = quatmultiply(W__dataS.Pelvis.ori, quatconj(xsensCalibSB.Pelvis.ori));
        qLankleEst0 = quatmultiply(W__dataS.L_LowLeg.ori, quatconj(xsensCalibSB.L_LowLeg.ori));
        qRankleEst0 = quatmultiply(W__dataS.R_LowLeg.ori, quatconj(xsensCalibSB.R_LowLeg.ori));
        qOri.w__sx.PELV = qPelvisEst0(sIdx:eIdx, :);
        qOri.w__sx.LTIB = qLankleEst0(sIdx:eIdx, :);
        qOri.w__sx.RTIB = qRankleEst0(sIdx:eIdx, :);
        
        % orientation of body from xsens
        qOri.w__x.PELV = W__dataX.qHips(sIdx:eIdx, :);
        qOri.w__x.LTIB = W__dataX.qLeftLeg(sIdx:eIdx, :);
        qOri.w__x.RTIB = W__dataX.qRightLeg(sIdx:eIdx, :);
       
        % Position, Velocity, Acceleration
        % gfrAcc from xsens
        W__xsensBody = W__dataX.togrBody(1:nSamples, {'name', 'xsens', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                             'xyzColor', {'m', 'y', 'c'}}); 
        vel = W__xsensBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
        acc = W__xsensBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});

        x0.w__x = [W__xsensBody.MIDPEL(sIdx,:) vel.MIDPEL(sIdx,:) zeros(1,4) ...
                   W__xsensBody.LTIO(sIdx,:) vel.LTIO(sIdx,:) zeros(1,4) ...
                   W__xsensBody.RTIO(sIdx,:) vel.RTIO(sIdx,:) zeros(1,4)]';  
        
        gfrAcc.w__x = {};
        gfrAcc.w__x.MP = acc.MIDPEL;
        gfrAcc.w__x.MP = gfrAcc.w__x.MP(sIdx:eIdx,:);
        gfrAcc.w__x.LA = acc.LTIO;
        gfrAcc.w__x.LA = gfrAcc.w__x.LA(sIdx:eIdx,:);
        gfrAcc.w__x.RA = acc.RTIO;
        gfrAcc.w__x.RA = gfrAcc.w__x.RA(sIdx:eIdx,:);
        
        % gfrAcc from sparse
        gfrAcc.w__sx = {};
        gfrAcc.w__sx.MP = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        gfrAcc.w__sx.MP = gfrAcc.w__sx.MP(sIdx:eIdx,:);
        gfrAcc.w__sx.LA = quatrotate(quatconj(W__dataS.L_LowLeg.ori), ...
                                W__dataS.L_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__sx.LA = gfrAcc.w__sx.LA(sIdx:eIdx,:);
        gfrAcc.w__sx.RA = quatrotate(quatconj(W__dataS.R_LowLeg.ori), ...
                                W__dataS.R_LowLeg.acc) - [0 0 9.81];
        gfrAcc.w__sx.RA = gfrAcc.w__sx.RA(sIdx:eIdx,:);
        
        % UWB measurements
        %  Simulate uwb measurement by generating pairwise combinations, using the
        %  origin of each bone segment as the root point
        umBuf = struct;
        umBuf.left_tibia_mid_pelvis = vecnorm((W__xsensBody.MIDPEL-W__xsensBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_mid_pelvis = umBuf.left_tibia_mid_pelvis(sIdx:eIdx,:);
        umBuf.mid_pelvis_right_tibia = vecnorm((W__xsensBody.MIDPEL-W__xsensBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.mid_pelvis_right_tibia = umBuf.mid_pelvis_right_tibia(sIdx:eIdx,:);
        umBuf.left_tibia_right_tibia = vecnorm((W__xsensBody.RTIO-W__xsensBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.left_tibia_right_tibia = umBuf.left_tibia_right_tibia(sIdx:eIdx,:);
        umBuf.LLeg = vecnorm((W__xsensBody.LFEP-W__xsensBody.LTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.LLeg = umBuf.LLeg(sIdx:eIdx,:);
        umBuf.RLeg = vecnorm((W__xsensBody.RFEP-W__xsensBody.RTIO), 2, 2) ...
            + normrnd(0, uwbMeasSigma, [nSamples, 1]);
        umBuf.RLeg = umBuf.RLeg(sIdx:eIdx,:);
        uwbMeas.w__x = umBuf;
        
        % debug purposes
        W__xsensBody = W__dataX.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__xsensBody = W__xsensBody.changeRefFrame('MIDPEL');
    end
    
    
    %% Save processing
    if ~strcmp(savedir, '')
        if ~isempty(dataX)
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', 'V__viconBody', 'W__xsensBody', ...
                 'gfrAcc', 'qOri', 'x0', 'uwbMeas', 'allIdx')
        else
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', 'V__viconBody', ...
                 'gfrAcc', 'qOri', 'x0', 'uwbMeas', 'allIdx')
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
        elseif ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
           strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
            csGfrAcc = gfrAcc.(strcat(cs.accData, cs.initSrc(end)));
        else
            csGfrAcc = gfrAcc.(cs.accData);
        end
        
        if strcmp(cs.oriData, 'w__s') || strcmp(cs.oriData, 'v__s')
            csQOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
            csWOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
        else
            csQOri = qOri.(cs.oriData);
            csWOri = qOri.(cs.oriData);
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
        elseif strcmp(cs.initSrc, 'v__v')
            csActBody = V__viconBody;
            csActBodyRel = PV__viconBody;
            d_pelvis = norm(V__dataV.RFEP(sIdx,:) - V__dataV.LFEP(sIdx,:));
            d_rfemur = norm(V__dataV.RFEP(sIdx,:) - V__dataV.RFEO(sIdx,:));
            d_lfemur = norm(V__dataV.LFEP(sIdx,:) - V__dataV.LFEO(sIdx,:));
            d_rtibia = norm(V__dataV.RFEO(sIdx,:) - V__dataV.RTIO(sIdx,:));
            d_ltibia = norm(V__dataV.LFEO(sIdx,:) - V__dataV.LTIO(sIdx,:));
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
        if strcmp(cs.stepDetection, 'av01')
            VAR_WIN  = floor(fs*cs.stepDetectWindow); % NUM_SAMPLES
            ACC_VAR_THRESH = cs.stepDetectThreshold;

            movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
            bIsStatMP = movVarAcc_pelvis < 0;
            movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
            bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
            movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
            bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
        elseif strcmp(cs.stepDetection, 'av02')
            if cs.accData(1) == 'w', csGfrAcc2 = gfrAcc.w__v;
            else, csGfrAcc2 = gfrAcc.v__v; end
            
            VAR_WIN  = floor(fs*cs.stepDetectWindow); % NUM_SAMPLES
            ACC_VAR_THRESH = cs.stepDetectThreshold;

            movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc2.MP .^2, 2)), VAR_WIN);
            bIsStatMP = movVarAcc_pelvis < 0;
            movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc2.LA .^2, 2)), VAR_WIN);
            bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
            movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc2.RA .^2, 2)), VAR_WIN);
            bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
        elseif strcmp(cs.stepDetection, 'av03')
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
            if strcmp(cs.est, 'ekfv3')

                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigmaQAccMP', cs.sigmaQAcc, ...
                    'sigmaQAccLA', cs.sigmaQAcc, 'sigmaQAccRA', cs.sigmaQAcc, ...
                    'sigmaUwbLLeg', cs.sigmaUwbLeg, 'sigmaUwbRLeg', cs.sigmaUwbLeg, ...
                    'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);
%                 disp(sprintf('%.2f %.2f', rad2deg(alphaLKmin), rad2deg(alphaRKmin)));
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.kf_3_kmus_v3( ...
                    csx0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwbMeas.(cs.initSrc), v3Options);

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
            elseif strcmp(cs.est, 'mpfv1')
                v3Options = struct('fs', fs, 'applyPred', cs.applyPred, ...
                                   'applyMeas', cs.applyMeas );
                lkangle = csActBody.calcJointAnglesLKnee(1);
                rkangle = csActBody.calcJointAnglesRKnee(1);
                csx02 = [csx0(1:6, 1); csx0(11:16, 1); csx0(21:26, 1); ...
                         csActBody.qRPV(1,:)'; csActBody.qLSK(1,:)'; ...
                         csActBody.qRSK(1,:)'; lkangle(2); rkangle(2)];
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.mpf_3_kmus_v1( ...
                    csx02, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwbMeas.(cs.initSrc), v3Options);
                
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
                   'LTIO', x_pos_v2(idx1, 7:9), ...
                   'RFEP', t_dat_v2.RFEP(idx1, :), ...
                   'RFEO', t_dat_v2.RFEO(idx1, :), ...
                   'RTIO', x_pos_v2(idx1, 13:15), ...
                   'qRPV', x_pos_v2(idx1, 19:22), ...
                   'qLTH', t_dat_v2.qLTH(idx1, :), ...
                   'qRTH', t_dat_v2.qRTH(idx1, :), ...
                   'qLSK', x_pos_v2(idx1, 23:26), ...
                   'qRSK', x_pos_v2(idx1, 27:30));
                estState = x_pos_v2;
                estState2 = t_dat_v2;
            elseif strcmp(cs.est, 'lieekfv1')
                v3Options = struct('fs', fs);
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.lieekf_3_kmus_v1( ...
                    csx0, cs.P, ...
                    csGfrAcc.MP, bIsStatMP, csQOri.PELV, csWOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, csWOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, csWOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwbMeas.(cs.initSrc), v3Options);

                idx1EndIdx = find(any(isnan(x_pos_v2.x), 2), 1);
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
                   'MIDPEL', squeeze(x_pos_v2.RPV(1:3,4,idx1))', ...
                   'LFEP', t_dat_v2.LFEP(idx1, :), ...
                   'LFEO', t_dat_v2.LFEO(idx1, :), ...
                   'LTIO', squeeze(x_pos_v2.LSK(1:3,4,idx1))', ...
                   'RFEP', t_dat_v2.RFEP(idx1, :), ...
                   'RFEO', t_dat_v2.RFEO(idx1, :), ...
                   'RTIO', squeeze(x_pos_v2.RSK(1:3,4,idx1))', ...
                   'qRPV', rotm2quat(x_pos_v2.RPV(1:3,1:3,idx1)), ...
                   'qLTH', t_dat_v2.qLTH(idx1, :), ...
                   'qRTH', t_dat_v2.qRTH(idx1, :), ...
                   'qLSK', rotm2quat(x_pos_v2.LSK(1:3,1:3,idx1)), ...
                   'qRSK', rotm2quat(x_pos_v2.RSK(1:3,1:3,idx1)));
               
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
            results0 = estBody2.diffRMSEandMean(csActBody2);
%         catch
%             runtime = cputime-t0;
%             results0 = csActBodyRel.diffRMSEandMean(nan);
%         end
        
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