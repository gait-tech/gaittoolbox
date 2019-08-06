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
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+Aw__sOw__sIw__v+Sav03+M76+C355"), ....
      struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P023+M111+C007"), ...
      struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P023+M111+C007"), ...
%       struct('file', 'S01-Trial-FigureofEight-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P001+M111+C007"), ...
%       struct('file', 'S01-Trial-Fivemin-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P001+M111+C007"), ...
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

addpath('btk');

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
    if strcmp(cs.oriData, 'w__s') || strcmp(cs.oriData, 'v__s')
        csQOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
        csBodyWOri = wbodyOri.(strcat(cs.oriData, cs.initSrc(end)));
    else
        csQOri = qOri.(cs.oriData);
        csBodyWOri = wbodyOri.(cs.oriData);
    end
        
    idx = allIdx.(cs.initSrc);
    
    eMarkers = struct();

    estBodyRel = estBody.changeRefFrame('MIDPEL');
    viconBodyRel = vb.changeRefFrame('MIDPEL');

    estBody2 = estBodyRel.toWorldFrame(vb.MIDPEL, estBody.qRPV);
    viconBody2 = viconBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
    
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
    vel = vb.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
    angvel = vb.calcSegAngVel({'qRPV', 'qLSK', 'qRSK'}, 'W');
    sensors.PELVVelRef = vel.MIDPEL;
    sensors.LANKVelRef = vel.LTIO;
    sensors.RANKVelRef = vel.RTIO;
    sensors.PELVAVelRef = angvel.qRPV;
    sensors.LANKAVelRef = angvel.qLSK;
    sensors.RANKAVelRef = angvel.qRSK;
    
    seq = 'YXZ';    
    [r2 r1 r3] = quat2angle(csQOri.PELV, seq); eul = rad2deg([r1 r2 r3]);
    sensors.PelvisAnglesXsens = eul;
    [r2 r1 r3] = quat2angle(csQOri.LTIB, seq); eul = rad2deg([r1 r2 r3]);
    sensors.LShankAnglesXsens = eul;
    [r2 r1 r3] = quat2angle(csQOri.RTIB, seq); eul = rad2deg([r1 r2 r3]);
    sensors.RShankAnglesXsens = eul;
%     segs = {'PELV', 'LTIB', 'RTIB'};
%     for i=1:length(segs)
%         n = segs{i};
%         n2 = sprintf('%sAVelXsens', n);
%         w = quatmultiply(quatconj(csQOri.(n)(1:end, :)), ...
%                                 csQOri.(n)([2:end end], :));
%         tmpIdx = w(:,1)<0;
%         w(tmpIdx,:) = -w(tmpIdx,:);
%         
% %         sensors.(n2) = 2*w(:,2:4)*fs;
%         sensors.(n2) = quatrotate(quatconj(csQOri.(n)), 2*w(:,2:4)*vb.fs);
%     end
    sensors.ePos = estBody2.calcDPos(viconBody2);
    sensors.eOri = estBody2.calcDOri(viconBody2);
    
    sensors = experiment.buildSensorStructFromDebug(sensors, estState, estState2, cs.est);
    
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
        csGfrAcc = gfrAcc.(eLabel);
        
        csNSamples = size(csGfrAcc.MP, 1);
        bIsStatMP = false(csNSamples, 1);
        bIsStatLA = revStepDetect.stepL(idx);
        bIsStatRA = revStepDetect.stepR(idx);
    end

    % velocity constraint debug
    sensors = experiment.buildVelCstrDebug(sensors, estBody);
    sensors = experiment.buildVelCstrDebug(sensors, vb, 'Ref');
    
    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, eMarkers);
%     estBody.exportc3d(sprintf('%s-EstOnly.c3d', targetname), sensors, ...
%                       false, bIsStatLA, bIsStatRA);
%     estBody.exportc3d(sprintf('%s03.c3d', targetname), sensors, ...
%                       vb, bIsStatLA, bIsStatRA, eMarkers, 3);
    
    % side by side
%     estBodyS2S = estBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
%     estBodyS2S = estBodyRel.toWorldFrame(estBody.MIDPEL, estBody.qRPV);
    viconBodyS2S = viconBodyRel.toWorldFrame(vb.MIDPEL+[0.5 0 0], vb.qRPV);
    estBody.exportc3d(sprintf('%s-SidebySide.c3d', targetname), sensors, viconBodyS2S, bIsStatLA, bIsStatRA, eMarkers);
    
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
%     estBodyRel.exportc3d(sprintf('%s-Rel.c3d', targetname), sensors, viconBodyRel, bIsStatLA, bIsStatRA, eMarkers);
%     estBodyRel.exportc3d(sprintf('%s-Rel02.c3d', targetname), sensors, viconBodyRel, bIsStatLA, bIsStatRA, eMarkers, 2);

    estBodyViconPelv = estBodyRel.toWorldFrame(vb.MIDPEL, vb.qRPV);
%     estBodyViconPelv.exportc3d(sprintf('%s-Vicon.c3d', targetname), sensors, vb, bIsStatLA, bIsStatRA, eMarkers);

%     vb.exportc3d(sprintf('%s-Vicon.c3d', targetname));
    
    %% Debug level c3d export
    [estBodyDebug, sensorsDebug] = experiment.buildgrBodyDebug(estBody, sensors, estState2, cs.est);
    n2 = estBodyDebug.nSamples;
    [viconBodyDebug, ~] = experiment.buildgrBodyDebug(vb, {}, false, 'vicon');

    bIsStatLADebug = repelem(bIsStatLA, 3, 1);
    bIsStatRADebug = repelem(bIsStatRA, 3, 1);

    estBodyDebugRel = estBodyDebug.changeRefFrame('MIDPEL');
    viconBodyDebugRel = viconBodyDebug.changeRefFrame('MIDPEL');

    estBodyDebug2 = estBodyDebugRel.toWorldFrame(viconBodyDebug.MIDPEL, estBodyDebug.qRPV);
    viconBodyDebug2 = viconBodyDebugRel.toWorldFrame(viconBodyDebug.MIDPEL, viconBodyDebug.qRPV);
    
    sensorsDebug.ePos = estBodyDebug2.calcDPos(viconBodyDebug2);
    sensorsDebug.eOri = estBodyDebug2.calcDOri(viconBodyDebug2);
    
    estBodyDebug.exportc3d(sprintf('%s-Debug.c3d', targetname), sensorsDebug, ...
                           viconBodyDebug, bIsStatLADebug, bIsStatRADebug);
%     estBodyDebugRel.exportc3d(sprintf('%s-RelDebug.c3d', targetname), sensorsDebug, ...
%                               viconBodyDebugRel, bIsStatLADebug, bIsStatRADebug);
end

function out = addAxis(q, p, qname)
    out = struct();
    d = 0.1;
    
    R = quat2rotm(q);
    out.(sprintf('%sX', qname)) = p+d*squeeze(R(:,1,:))';
    out.(sprintf('%sY', qname)) = p+d*squeeze(R(:,2,:))';
    out.(sprintf('%sZ', qname)) = p+d*squeeze(R(:,3,:))';
end