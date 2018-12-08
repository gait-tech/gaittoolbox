% motion list
list = {
    struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dxx+M02+C178'), ...
    struct('subj', 's1', 'act', 'walking3', 'algo', 'TC1+Dxx+M02+C178'), ...
%     struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dxx+M02+C202'), ...
%     struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dxx+M02+C172'), ...
%     struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dvx+M02+C202'), ...
%     struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dxv+M02+C202'), ...
%     struct('subj', 's1', 'act', 'walking1', 'algo', 'TC1+Dvv+M02+C202'), ...
};

for i=1:length(list)
    dataSfname = sprintf('totalcapture/gyroMag/%s/%s_Xsens_AuxFields.sensors', ...
                         list{i}.subj, list{i}.act);
    load(sprintf('totalcapture/explore/tcd-%s-%s-debug.mat', ...
                 list{i}.subj, list{i}.act));
    load(sprintf('totalcapture/explore/tcd-%s-%s-%s.mat', ...
                 list{i}.subj, list{i}.act, list{i}.algo));
    targetname = sprintf('explore_output/tcd-%s-%s-%s', ...
                         list{i}.subj, list{i}.act, list{i}.algo);
    
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
    viconBodyRel = actBody.changeRefFrame('MIDPEL');

    d = norm(estBody.MIDPEL(1,:) - estBody.LFEP(1,:))*0.5;

    dataS = mocapdb.XsensBody.loadSensorFile(dataSfname);
    dataS = dataS.getSubset(idx);

    % sensors = struct();
    sensors = dataS.exportRawMeasurementAsStruct({'Pelvis', 'L_LowLeg', 'R_LowLeg'}, ...
                    {'PELV', 'LANK', 'RANK'});
    sensors.PELVFreeAccX = gfrAcc.x.MP(:,1);
    sensors.PELVFreeAccY = gfrAcc.x.MP(:,2);
    sensors.PELVFreeAccZ = gfrAcc.x.MP(:,3);
    sensors.LANKFreeAccX = gfrAcc.x.LA(:,1);
    sensors.LANKFreeAccY = gfrAcc.x.LA(:,2);
    sensors.LANKFreeAccZ = gfrAcc.x.LA(:,3);
    sensors.RANKFreeAccX = gfrAcc.x.RA(:,1);
    sensors.RANKFreeAccY = gfrAcc.x.RA(:,2);
    sensors.RANKFreeAccZ = gfrAcc.x.RA(:,3);

    sensors.PELVFreeAccRefX = gfrAcc.v.MP(:,1);
    sensors.PELVFreeAccRefY = gfrAcc.v.MP(:,2);
    sensors.PELVFreeAccRefZ = gfrAcc.v.MP(:,3);
    sensors.LANKFreeAccRefX = gfrAcc.v.LA(:,1);
    sensors.LANKFreeAccRefY = gfrAcc.v.LA(:,2);
    sensors.LANKFreeAccRefZ = gfrAcc.v.LA(:,3);
    sensors.RANKFreeAccRefX = gfrAcc.v.RA(:,1);
    sensors.RANKFreeAccRefY = gfrAcc.v.RA(:,2);
    sensors.RANKFreeAccRefZ = gfrAcc.v.RA(:,3);

    % step detection
    fs = estBody.fs;
    VAR_WIN  = floor(fs*0.25); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;
    csGfrAcc = gfrAcc.x;

    movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
    bIsStatMP = movVarAcc_pelvis < 0;
    movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
    bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
    movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
    bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;

    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, actBody, bIsStatLA, bIsStatRA, eMarkers);
    
    % side by side
    estBodyS2S = estBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
    viconBodyS2S = viconBodyRel.toWorldFrame(actBody.MIDPEL+[0.5 0 0], actBody.qRPV);
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

    estBodyViconPelv = estBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
    estBodyViconPelv.exportc3d(sprintf('%s-Vicon.c3d', targetname), sensors, actBody, bIsStatLA, bIsStatRA, eMarkers);

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

    viconBodyDebug = actBody.copy();
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
    % viconBody.exportc3d(sprintf('%svicon.c3d', targetname));
end

function out = addAxis(q, p, qname)
    out = struct();
    d = 0.1;
    
    R = quat2rotm(q);
    out.(sprintf('%sX', qname)) = p+d*squeeze(R(:,1,:))';
    out.(sprintf('%sY', qname)) = p+d*squeeze(R(:,2,:))';
    out.(sprintf('%sZ', qname)) = p+d*squeeze(R(:,3,:))';
end