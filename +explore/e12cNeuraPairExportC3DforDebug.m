% motion list
list = {
    % vicon result investigation
      struct('file', 'S01-Trial-HighKneeJog-1', ...
             'algo1', "NS2+Aw__sOw__sIw__v+Sav03+M302+C001", ...
             'algo2', "NS2+Aw__sOw__sIw__v+Sav03+M302+C351", ...
             'name', "NS2+Aw__sOw__sIw__v+Sav03+M302+C001vs351"), ...
%       struct('file', 'S02-Trial-Walk-2', ...
%              'algo1', "NS1+Aw__sOw__sIw__v+Sav01+M76+C355", ...
%              'algo2', "NS1+Aw__sOw__sIw__v+Sav01+M76+C375", ...
%              'name', "NS1+Aw__sOw__sIw__v+Sav01+M76+C355vs375"), ...
%       struct('file', 'S02-Trial-Walk-2', ...
%              'algo1', "NS1+Aw__sOw__sIw__v+Sav01+M76+C355", ...
%              'algo2', "NS1+Aw__sOw__sIw__v+Sav01+M76+C353", ...
%              'name', "NS1+Aw__sOw__sIw__v+Sav01+M76+C355vs353"), ...
};
ns = "NS2";
for lIdx=1:length(list)
    dataSfname = sprintf('neura-sparse01/imu/%s', list{lIdx}.file);
    load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, list{lIdx}.file));
    
    targetname = sprintf('explore_output/%s-%s-%s', ns, list{lIdx}.file, list{lIdx}.name);
    
    %% algorithm 1
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo1));
    
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
    estBody1 = estBody.copy();
    estBody1Rel = estBody1.changeRefFrame('MIDPEL');
    viconBodyRel = vb.changeRefFrame('MIDPEL');

    d = norm(estBody1.MIDPEL(1,:) - estBody1.LFEP(1,:))*0.5;

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
    fs = estBody1.fs;
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
    end
   
    % side by side
%     estBody1S2S = estBody1Rel.toWorldFrame(vb.MIDPEL, vb.qRPV);
    estBody1S2S = estBody1Rel.toWorldFrame(vb.MIDPEL, estBody1.qRPV);
    viconBodyS2S = viconBodyRel.toWorldFrame(vb.MIDPEL+[0.5 0 0], vb.qRPV);
    
    estBody1Debug = estBody1.copy();
    n2 = estBody1.nSamples*3;
    estBody1Debug.MIDPEL = repelem(estBody1Debug.MIDPEL, 3, 1);
    estBody1Debug.MIDPEL(1:3:n2, :) = estState2.predState(:,1:3);
    estBody1Debug.MIDPEL(2:3:n2, :) = estState2.zuptState(:,1:3);
    estBody1Debug.LTIO = repelem(estBody1Debug.LTIO, 3, 1);
    estBody1Debug.LTIO(1:3:n2, :) = estState2.predState(:,11:13);
    estBody1Debug.LTIO(2:3:n2, :) = estState2.zuptState(:,11:13);
    estBody1Debug.RTIO = repelem(estBody1Debug.RTIO, 3, 1);
    estBody1Debug.RTIO(1:3:n2, :) = estState2.predState(:,21:23);
    estBody1Debug.RTIO(2:3:n2, :) = estState2.zuptState(:,21:23);

    estBody1Debug.qRPV = repelem(estBody1Debug.qRPV, 3, 1);
    estBody1Debug.qLSK = repelem(estBody1Debug.qLSK, 3, 1);
    estBody1Debug.qRSK = repelem(estBody1Debug.qRSK, 3, 1);

    v = quat2rotm(estBody1Debug.qLSK); v = squeeze(v(:,3,:))';
    estBody1Debug.LFEO = estBody1Debug.LTIO + estBody1Debug.getLShankLength()*v;
    v = quat2rotm(estBody1Debug.qRSK); v = squeeze(v(:,3,:))';
    estBody1Debug.RFEO = estBody1Debug.RTIO + estBody1Debug.getRShankLength()*v;
    v = quat2rotm(estBody1Debug.qRPV); v = squeeze(v(:,2,:))';
    estBody1Debug.LFEP = estBody1Debug.MIDPEL + estBody1Debug.getPelvisLength()/2*v;
    estBody1Debug.RFEP = estBody1Debug.MIDPEL - estBody1Debug.getPelvisLength()/2*v;

    v = zeros(3, 3, n2);
    z = (estBody1Debug.LFEP-estBody1Debug.LFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBody1Debug.qLSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBody1Debug.qLTH = rotm2quat(v);

    v = zeros(3, 3, n2);
    z = (estBody1Debug.RFEP-estBody1Debug.RFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBody1Debug.qRSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBody1Debug.qRTH = rotm2quat(v);

    estBody1Debug.nSamples = estBody1Debug.nSamples*3;
    estBody1Debug.fs = estBody1Debug.fs*3;

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
    
    %% algorithm 2
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo2));
    estBody2 = estBody.copy();
    estBody2Rel = estBody2.changeRefFrame('MIDPEL');
    estBody2S2S = estBody2Rel.toWorldFrame(vb.MIDPEL+[0.5 0 0], estBody2.qRPV);
    estBody2Debug = estBody2.copy();
    n2 = estBody1.nSamples*3;
    estBody2Debug.MIDPEL = repelem(estBody2Debug.MIDPEL, 3, 1);
    estBody2Debug.MIDPEL(1:3:n2, :) = estState2.predState(:,1:3);
    estBody2Debug.MIDPEL(2:3:n2, :) = estState2.zuptState(:,1:3);
    estBody2Debug.LTIO = repelem(estBody2Debug.LTIO, 3, 1);
    estBody2Debug.LTIO(1:3:n2, :) = estState2.predState(:,11:13);
    estBody2Debug.LTIO(2:3:n2, :) = estState2.zuptState(:,11:13);
    estBody2Debug.RTIO = repelem(estBody2Debug.RTIO, 3, 1);
    estBody2Debug.RTIO(1:3:n2, :) = estState2.predState(:,21:23);
    estBody2Debug.RTIO(2:3:n2, :) = estState2.zuptState(:,21:23);

    estBody2Debug.qRPV = repelem(estBody2Debug.qRPV, 3, 1);
    estBody2Debug.qLSK = repelem(estBody2Debug.qLSK, 3, 1);
    estBody2Debug.qRSK = repelem(estBody2Debug.qRSK, 3, 1);

    v = quat2rotm(estBody2Debug.qLSK); v = squeeze(v(:,3,:))';
    estBody2Debug.LFEO = estBody2Debug.LTIO + estBody2Debug.getLShankLength()*v;
    v = quat2rotm(estBody2Debug.qRSK); v = squeeze(v(:,3,:))';
    estBody2Debug.RFEO = estBody2Debug.RTIO + estBody2Debug.getRShankLength()*v;
    v = quat2rotm(estBody2Debug.qRPV); v = squeeze(v(:,2,:))';
    estBody2Debug.LFEP = estBody2Debug.MIDPEL + estBody2Debug.getPelvisLength()/2*v;
    estBody2Debug.RFEP = estBody2Debug.MIDPEL - estBody2Debug.getPelvisLength()/2*v;

    v = zeros(3, 3, n2);
    z = (estBody2Debug.LFEP-estBody2Debug.LFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBody2Debug.qLSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBody2Debug.qLTH = rotm2quat(v);

    v = zeros(3, 3, n2);
    z = (estBody2Debug.RFEP-estBody2Debug.RFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(estBody2Debug.qRSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    estBody2Debug.qRTH = rotm2quat(v);

    estBody2Debug.nSamples = estBody2Debug.nSamples*3;
    estBody2Debug.fs = estBody2Debug.fs*3;
    
    sensors.PELVVelRef = estState(:, 4:6);
    sensors.LANKVelRef = estState(:, 14:16);
    sensors.RANKVelRef = estState(:, 24:26);
    
    fn = {"PELVVelRef", "LANKVelRef", "RANKVelRef"};
    for i=1:length(fn)
        n = fn{i};
        sensorsDebug.(n) = repelem(sensors.(n), 3, 1);
    end
    sensorsDebug.PELVVelRef(1:3:n2, :) = estState2.predState(:,  4: 6);
    sensorsDebug.PELVVelRef(2:3:n2, :) = estState2.zuptState(:,  4: 6);
    sensorsDebug.LANKVelRef(1:3:n2, :) = estState2.predState(:, 14:16);
    sensorsDebug.LANKVelRef(2:3:n2, :) = estState2.zuptState(:, 14:16);
    sensorsDebug.RANKVelRef(1:3:n2, :) = estState2.predState(:, 24:26);
    sensorsDebug.RANKVelRef(2:3:n2, :) = estState2.zuptState(:, 24:26);
    
    %% export
    estBody1.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                       estBody2, bIsStatLA, bIsStatRA, eMarkers);
    estBody1S2S.exportc3d(sprintf('%s-SidebySide.c3d', targetname), ...
                sensors, estBody2S2S, bIsStatLA, bIsStatRA, eMarkers);
    estBody1Debug.exportc3d(sprintf('%s-Debug.c3d', targetname), sensorsDebug, ...
                           estBody2Debug, bIsStatLADebug, bIsStatRADebug);
end

function out = addAxis(q, p, qname)
    out = struct();
    d = 0.1;
    
    R = quat2rotm(q);
    out.(sprintf('%sX', qname)) = p+d*squeeze(R(:,1,:))';
    out.(sprintf('%sY', qname)) = p+d*squeeze(R(:,2,:))';
    out.(sprintf('%sZ', qname)) = p+d*squeeze(R(:,3,:))';
end