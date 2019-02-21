% motion list
list = {
    % vicon result investigation
%     struct('file', 'S10-Trial-Static-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S02-Trial-Static-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S06-Trial-Static-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S05-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S08-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S02-Trial-JumpingJacks-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S04-Trial-JumpingJacks-2', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
%     struct('file', 'S01-Trial-JumpingJacks-1', 'algo', 'NS1+Av__sOv__sIv__v+Sav01+M02+C176'), ...
    % xsens result investigation
%     struct('file', 'S10-Trial-Static-1', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S02-Trial-Static-1', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S03-Trial-Static-1', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S02-Trial-Walk-2', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S08-Trial-Walk-2', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S06-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S09-Trial-JumpingJacks-2', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S04-Trial-JumpingJacks-2', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%     struct('file', 'S01-Trial-JumpingJacks-1', 'algo', 'NS1+Aw__sOw__sIw__x+Sav01+M02+C176'), ...
%      % debug
%       struct('file', 'S02-Trial-Walk-2', 'algo', "NS1+Aw__sOw__sIw__v+Sav01+M70+C355"), ...
%       struct('file', 'S03-Trial-Walk-1', 'algo', "NS1+Aw__sOw__sIw__v+Sav01+M76+C455"), ...
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+Aw__sOw__sIw__v+Sav01+M76+C455"), ...
      struct('file', 'S07-Trial-Walk-1', 'algo', "NS2+Aw__sOw__sIw__v+Sav01+M76+C355"), ...
      struct('file', 'S07-Trial-Walk-1', 'algo', "NS2+Aw__sOw__sIw__v+Sav03+M76+C355"), ...
%       struct('file', 'S02-Trial-Walk-2', 'algo', "NS1+Aw__sOw__sIw__v+Sav01+M72+C355"), ...
%       struct('file', 'S02-Trial-Walk-2', 'algo', "NS1+Aw__sOw__sIw__v+Sav01+M74+C355"), ...
%       struct('file', 'S02-Trial-Walk-2', 'algo', "NS1+Aw__sOw__sIw__v+Sav01+M76+C355"), ...
    % london demo
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M00+C000"), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M02+C000"), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M02+C273"), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M02+C205"), ...
    % uwb test
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M102+C201"), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M102+C202"), ...
%     struct('file', 'S07-Trial-Walk-1', 'algo', "NS1+Av__sOv__sIv__v+Sav01+M102+C175"), ...
};

for lIdx=1:length(list)
    dataSfname = sprintf('neura-sparse01/imu/%s', list{lIdx}.file);
    ns = extractBetween(list{lIdx}.algo, 1, 3);
    
    load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, list{lIdx}.file));
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo));
    targetname = sprintf('explore_output/%s-%s-%s', ns, list{lIdx}.file, list{lIdx}.algo);
    
    if strcmp(cs.stepDetection, 'av03')
        revStepDetect = readtable(sprintf('neura-sparse01/step-detect/%s-revStepDetect.csv', list{lIdx}.file));
    end
    
    if cs.initSrc == 'w__v'
        aLabel = 'w__v';
        vb = W__viconBody;
    elseif cs.initSrc == 'v__v'
        aLabel = 'v__v';
        vb = V__viconBody;
    else
        aLabel = 'w__x';
        vb = W__xsensBody;
    end
    if ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
       strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
        eLabel = strcat(cs.accData, cs.initSrc(end));
    else
        eLabel = cs.accData;
    end
    idx = allIdx.(cs.initSrc);
    
    eMarkers = struct();

    % estBodyPred = estBody.copy();
    % estBodyPred.MIDPEL = estState2.predState(:,1:3);
    % estBodyPred.LTIO = estState2.predState(:,11:13);
    % estBodyPred.RTIO = estState2.predState(:,21:23);
    % estBodyPredRel = estBodyPred.changeRefFrame('MIDPEL');
    % 
    % eMarkers.MIDPELPred = estBodyPredRel.MIDPEL;
    % eMarkers.LTIOPred = estBodyPredRel.LTIO;
    % eMarkers.RTIOPred = estBodyPredRel.RTIO;
    % 
    % estBodyMeas = estBody.copy();
    % estBodyMeas.MIDPEL = estState2.zuptState(:,1:3);
    % estBodyMeas.LTIO = estState2.zuptState(:,11:13);
    % estBodyMeas.RTIO = estState2.zuptState(:,21:23);
    % estBodyMeasRel = estBodyMeas.changeRefFrame('MIDPEL');
    % 
    % eMarkers.MIDPELMeas = estBodyMeasRel.MIDPEL;
    % eMarkers.LTIOMeas = estBodyMeasRel.LTIO;
    % eMarkers.RTIOMeas = estBodyMeasRel.RTIO;

    estBodyRel = estBody.changeRefFrame('MIDPEL');
    viconBodyRel = vb.changeRefFrame('MIDPEL');

    d = norm(estBody.MIDPEL(1,:) - estBody.LFEP(1,:))*0.5;

    options = struct('Pelvis', '00B40B91', ...
        'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
        'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
        'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
    dataS = mocapdb.XsensBody.loadMTExport(dataSfname, options);
    dataS = dataS.getSubset(idx);

%     sensors = struct();
    sensors = dataS.exportRawMeasurementAsStruct({'Pelvis', 'L_LowLeg', 'R_LowLeg'}, ...
                    {'PELV', 'LANK', 'RANK'});
    sensors.PELVFreeAcc = gfrAcc.(eLabel).MP;
    sensors.LANKFreeAcc = gfrAcc.(eLabel).LA;
    sensors.RANKFreeAcc = gfrAcc.(eLabel).RA;
    sensors.PELVFreeAccRef = gfrAcc.(aLabel).MP;
    sensors.LANKFreeAccRef = gfrAcc.(aLabel).LA;
    sensors.RANKFreeAccRef = gfrAcc.(aLabel).RA;
    sensors.PELVVel = estState(:, 4:6);
    sensors.LANKVel = estState(:, 14:16);
    sensors.RANKVel = estState(:, 24:26);
    vel = vb.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
    sensors.PELVVelRef = vel.MIDPEL;
    sensors.LANKVelRef = vel.LTIO;
    sensors.RANKVelRef = vel.RTIO;
    
    % step detection
    fs = estBody.fs;
    VAR_WIN  = floor(fs*0.25); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;

    if strcmp(cs.stepDetection, 'av01')
        csGfrAcc = gfrAcc.(eLabel);
        
        movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
        bIsStatMP = movVarAcc_pelvis < 0;
        movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
        bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
        movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
        bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
    elseif strcmp(cs.stepDetection, 'av02')
        if cs.accData(1) == 'w', csGfrAcc2 = gfrAcc.w__v;
        else, csGfrAcc2 = gfrAcc.v__v; end

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
    end

    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, eMarkers);
    estBody.exportc3d(sprintf('%s-EstOnly.c3d', targetname), sensors, ...
                      false, bIsStatLA, bIsStatRA);
    estBody.exportc3d(sprintf('%s03.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, eMarkers, 3);
    
    % side by side
%     estBodyS2S = estBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
    estBodyS2S = estBodyRel.toWorldFrame(vb.MIDPEL, estBody.qRPV);
    viconBodyS2S = viconBodyRel.toWorldFrame(vb.MIDPEL+[0.5 0 0], vb.qRPV);
    estBodyS2S.exportc3d(sprintf('%s-SidebySide.c3d', targetname), sensors, viconBodyS2S, bIsStatLA, bIsStatRA, eMarkers);
    
    tmp = addAxis(quatmultiply(quatconj(estBody.qRPV), dataS.Pelvis.ori), ...
                  estBodyRel.MIDPEL, 'qRPVSens');
    f = fieldnames(tmp);
    for i = 1:length(f)
        eMarkers.(f{i}) = tmp.(f{i});
    end
    tmp = addAxis(quatmultiply(quatconj(estBody.qRPV), dataS.L_LowLeg.ori), ...
                  estBodyRel.LTIO, 'qLSKSens');
    f = fieldnames(tmp);
    for i = 1:length(f)
        eMarkers.(f{i}) = tmp.(f{i});
    end
    tmp = addAxis(quatmultiply(quatconj(estBody.qRPV), dataS.R_LowLeg.ori), ...
                  estBodyRel.RTIO, 'qRSKSens');
    f = fieldnames(tmp);
    for i = 1:length(f)
        eMarkers.(f{i}) = tmp.(f{i});
    end
    estBodyRel.exportc3d(sprintf('%s-Rel.c3d', targetname), sensors, viconBodyRel, bIsStatLA, bIsStatRA, eMarkers);
    estBodyRel.exportc3d(sprintf('%s-Rel02.c3d', targetname), sensors, viconBodyRel, bIsStatLA, bIsStatRA, eMarkers, 2);

    estBodyViconPelv = estBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
    estBodyViconPelv.exportc3d(sprintf('%s-Vicon.c3d', targetname), sensors, vb, bIsStatLA, bIsStatRA, eMarkers);

    estBodyDebug = estBody.copy();
    n2 = estBody.nSamples*3;
    estBodyDebug.MIDPEL = repelem(estBodyDebug.MIDPEL, 3, 1);
    estBodyDebug.MIDPEL(1:3:n2, :) = estState2.predState(:,1:3);
    estBodyDebug.MIDPEL(2:3:n2, :) = estState2.zuptState(:,1:3);
    estBodyDebug.LTIO = repelem(estBodyDebug.LTIO, 3, 1);
    estBodyDebug.LTIO(1:3:n2, :) = estState2.predState(:,11:13);
    estBodyDebug.LTIO(2:3:n2, :) = estState2.zuptState(:,11:13);
    estBodyDebug.RTIO = repelem(estBodyDebug.RTIO, 3, 1);
    estBodyDebug.RTIO(1:3:n2, :) = estState2.predState(:,21:23);
    estBodyDebug.RTIO(2:3:n2, :) = estState2.zuptState(:,21:23);

    estBodyDebug.qRPV = repelem(estBodyDebug.qRPV, 3, 1);
    estBodyDebug.qLSK = repelem(estBodyDebug.qLSK, 3, 1);
    estBodyDebug.qRSK = repelem(estBodyDebug.qRSK, 3, 1);

    v = quat2rotm(estBodyDebug.qLSK); v = squeeze(v(:,3,:))';
    estBodyDebug.LFEO = estBodyDebug.LTIO + estBodyDebug.getLShankLength()*v;
    v = quat2rotm(estBodyDebug.qRSK); v = squeeze(v(:,3,:))';
    estBodyDebug.RFEO = estBodyDebug.RTIO + estBodyDebug.getRShankLength()*v;
    v = quat2rotm(estBodyDebug.qRPV); v = squeeze(v(:,2,:))';
    estBodyDebug.LFEP = estBodyDebug.MIDPEL + estBodyDebug.getPelvisLength()/2*v;
    estBodyDebug.RFEP = estBodyDebug.MIDPEL - estBodyDebug.getPelvisLength()/2*v;

    v = zeros(3, 3, n2);
    z = (estBodyDebug.LFEP-estBodyDebug.LFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBodyDebug.qLSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBodyDebug.qLTH = rotm2quat(v);

    v = zeros(3, 3, n2);
    z = (estBodyDebug.RFEP-estBodyDebug.RFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBodyDebug.qRSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBodyDebug.qRTH = rotm2quat(v);

    estBodyDebug.nSamples = estBodyDebug.nSamples*3;
    estBodyDebug.fs = estBodyDebug.fs*3;

    viconBodyDebug = vb.copy();
    for i=1:length(viconBodyDebug.posList)
        n = viconBodyDebug.posList{i};
        viconBodyDebug.(n) = repelem(viconBodyDebug.(n), 3, 1);
    end
    for i=1:length(viconBodyDebug.oriList)
        n = viconBodyDebug.oriList{i};
        viconBodyDebug.(n) = repelem(viconBodyDebug.(n), 3, 1);
    end
    viconBodyDebug.nSamples = viconBodyDebug.nSamples*3;

    bIsStatLADebug = repelem(bIsStatLA, 3, 1);
    bIsStatRADebug = repelem(bIsStatRA, 3, 1);

    estBodyDebugRel = estBodyDebug.changeRefFrame('MIDPEL');
    viconBodyDebugRel = viconBodyDebug.changeRefFrame('MIDPEL');

    sensorsDebug = struct();
    fn = fieldnames(sensors);
    for i=1:length(fn)
        n = fn{i};
        sensorsDebug.(n) = repelem(sensors.(n), 3, 1);
    end
    sensorsDebug.PELVVel(1:3:n2, :) = estState2.predState(:,  4: 6);
    sensorsDebug.PELVVel(2:3:n2, :) = estState2.zuptState(:,  4: 6);
    sensorsDebug.LANKVel(1:3:n2, :) = estState2.predState(:, 14:16);
    sensorsDebug.LANKVel(2:3:n2, :) = estState2.zuptState(:, 14:16);
    sensorsDebug.RANKVel(1:3:n2, :) = estState2.predState(:, 24:26);
    sensorsDebug.RANKVel(2:3:n2, :) = estState2.zuptState(:, 24:26);
    estBodyDebug.exportc3d(sprintf('%s-Debug.c3d', targetname), sensorsDebug, ...
                           viconBodyDebug, bIsStatLADebug, bIsStatRADebug);
    estBodyDebugRel.exportc3d(sprintf('%s-RelDebug.c3d', targetname), sensorsDebug, ...
                              viconBodyDebugRel, bIsStatLADebug, bIsStatRADebug);
    vb.exportc3d(sprintf('%s-Vicon.c3d', targetname));
end

function out = addAxis(q, p, qname)
    out = struct();
    d = 0.1;
    
    R = quat2rotm(q);
    out.(sprintf('%sX', qname)) = p+d*squeeze(R(:,1,:))';
    out.(sprintf('%sY', qname)) = p+d*squeeze(R(:,2,:))';
    out.(sprintf('%sZ', qname)) = p+d*squeeze(R(:,3,:))';
end