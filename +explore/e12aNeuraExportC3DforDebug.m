% motion list
dir = 'neura-sparse01';
dataList = readtable(sprintf('%s/data-list-v2.csv', dir));

list = {
%      % debug
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P021+M105+C007"), ....
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P021+M305+C007"), ...
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P021+M505+C007"), ...
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P021+M105+C017"), ...
%       struct('file', 'S01-Trial-Walk-1', 'algo', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P021+M105+C027"), ...
};

listN = size(list, 1);
for i = 135
    n = table2struct(dataList(i, :));
    name = sprintf("%s-%s", n.subj, n.act);
    for sI = [struct('algo', 'lieekfv1', 'pI', 1, 'mI', 125, 'cI', 7)] 
         list{listN+1} = struct('file', name, ...
             'algo', sprintf("NS2+%s+Aw__sOw__sIw__v+Sav03+P%03d+M%03d+C%03d", ...
                     sI.algo, sI.pI, sI.mI, sI.cI));
         listN = listN+1;
    end
end

addpath('btk');

targetdir = "explore_output";

for lIdx=1:length(list)
    dataSfname = sprintf('neura-sparse01/imu/%s', list{lIdx}.file);
    ns = extractBetween(list{lIdx}.algo, 1, 3);
    
    load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, list{lIdx}.file));
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo));
    targetname = sprintf("%s/%s-%s-%s", ...
                         targetdir, ns, list{lIdx}.file, list{lIdx}.algo);
    
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
    sensors = experiment.buildVelCstrDebug(sensors, estBody, estState, cs.est);
    sensors = experiment.buildVelCstrDebug(sensors, vb, estState, cs.est, 'Ref');
    
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
%     estBody.exportc3d(sprintf('%s-SidebySide.c3d', targetname), sensors, viconBodyS2S, bIsStatLA, bIsStatRA, eMarkers);
    
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
    
    estStateDebug = experiment.buildStateDebug(estState, estState2, cs.est);
    sensorsDebug = experiment.buildVelCstrDebug(sensorsDebug, estBodyDebug, estStateDebug, cs.est);
    
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