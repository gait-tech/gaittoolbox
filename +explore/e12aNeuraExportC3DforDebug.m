% motion list
list = {
%     struct('file', 'S01-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+M02+C201'), ...
    struct('file', 'S05-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+M02+C201'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+M02+C203'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Av__sOv__sIv__v+M02+C203'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Av__sOw__sIw__v+M02+C201'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Av__vOv__vIv__v+M02+C201'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C201'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C202'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C172'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C202'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C202'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'NS1+Aw__sOw__sIw__v+M02+C202'), ...
%     struct('file', 'S01-Trial-Walk-1', 'algo', 'Nssv+M02+C201'), ...
%     struct('file', 'S02-Trial-Walk-1', 'algo', 'Nssv+M02+C201'), ...
%     struct('file', 'S03-Trial-Walk-1', 'algo', 'Nssv+M02+C201'), ...
};

for lIdx=1:length(list)
    dataSfname = sprintf('neura-sparse01/imu/%s', list{lIdx}.file);
    load(sprintf('neura-sparse01/explore/neura-%s-debug.mat', list{lIdx}.file));
    load(sprintf('neura-sparse01/explore/neura-%s-%s.mat', list{lIdx}.file, list{lIdx}.algo));
    
    targetname = sprintf('explore_output/neura-%s-%s', list{lIdx}.file, list{lIdx}.algo);
    
    if cs.accData == 'w__s'
        eLabel = 'w__s';
        aLabel = 'w__v';
        vb = W__viconBody;
    else
        eLabel = 'v__s';
        aLabel = 'v__v';
        vb = V__viconBody;
    end
    
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
%     sensors.PELVFreeAccX = gfrAcc.(eLabel).MP(:,1);
%     sensors.PELVFreeAccY = gfrAcc.(eLabel).MP(:,2);
%     sensors.PELVFreeAccZ = gfrAcc.(eLabel).MP(:,3);
%     sensors.LANKFreeAccX = gfrAcc.(eLabel).LA(:,1);
%     sensors.LANKFreeAccY = gfrAcc.(eLabel).LA(:,2);
%     sensors.LANKFreeAccZ = gfrAcc.(eLabel).LA(:,3);
%     sensors.RANKFreeAccX = gfrAcc.(eLabel).RA(:,1);
%     sensors.RANKFreeAccY = gfrAcc.(eLabel).RA(:,2);
%     sensors.RANKFreeAccZ = gfrAcc.(eLabel).RA(:,3);
% 
%     sensors.PELVFreeAccRefX = gfrAcc.(aLabel).MP(:,1);
%     sensors.PELVFreeAccRefY = gfrAcc.(aLabel).MP(:,2);
%     sensors.PELVFreeAccRefZ = gfrAcc.(aLabel).MP(:,3);
%     sensors.LANKFreeAccRefX = gfrAcc.(aLabel).LA(:,1);
%     sensors.LANKFreeAccRefY = gfrAcc.(aLabel).LA(:,2);
%     sensors.LANKFreeAccRefZ = gfrAcc.(aLabel).LA(:,3);
%     sensors.RANKFreeAccRefX = gfrAcc.(aLabel).RA(:,1);
%     sensors.RANKFreeAccRefY = gfrAcc.(aLabel).RA(:,2);
%     sensors.RANKFreeAccRefZ = gfrAcc.(aLabel).RA(:,3);
    sensors.PELVFreeAcc = gfrAcc.(eLabel).MP;
    sensors.LANKFreeAcc = gfrAcc.(eLabel).LA;
    sensors.RANKFreeAcc = gfrAcc.(eLabel).RA;
    sensors.PELVFreeAccRef = gfrAcc.(aLabel).MP;
    sensors.LANKFreeAccRef = gfrAcc.(aLabel).LA;
    sensors.RANKFreeAccRef = gfrAcc.(aLabel).RA;
    
    % step detection
    fs = estBody.fs;
    VAR_WIN  = floor(fs*0.25); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;
    csGfrAcc = gfrAcc.(eLabel);

    movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
    bIsStatMP = movVarAcc_pelvis < 0;
    movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
    bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
    bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;

    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, eMarkers);
    estBody.exportc3d(sprintf('%s03.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, eMarkers, 3);
    
    % side by side
    estBodyS2S = estBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
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

    viconBodyDebug = W__viconBody.copy();
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