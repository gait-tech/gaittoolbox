% motion list
list = {
%       struct('file', 'S01-Trial-Walk-1', ...
%              'algo1', "NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P001+M104+C000", ...
%              'algo2', "NS2+lieekfv2+Aw__sOw__sIw__v+Sav03+P001+M104+C000", ...
%              'name', "NS2+lgekfv1xv2+Aw__sOw__sIw__v+Sav03+P001+M104+C000"), ...
      
};
for sI = [struct('pI', 1, 'mI', 135, 'cI', 7), ...
              ] 
    list{end+1} = struct('file', 'S01-Trial-Walk-1', ...
         'algo1', sprintf("NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P%03d+M%03d+C%03d", ...
                          sI.pI, sI.mI, sI.cI), ...
         'algo2', sprintf("NS2+lieekfv1+Aw__sOw__sIw__v+Sav03+P%03d+M%03d+C%03d", ...
                          sI.pI, sI.mI-10, sI.cI), ...
         'name', sprintf("NS2+lieekfv1xv2+Aw__sOw__sIw__v+Sav03+P%03d+M%03d+C%03d", ...
                         sI.pI, sI.mI, sI.cI) );
end
ns = "NS2";
% targetdir = "C:/Users/z5151460/OneDrive - UNSW/Thesis - Sparse Mocap/Goal02-LGCEKF/Analysis-20190819-SE3xR6v2VSSO3xR9v1/c3d";
targetdir = "explore_output";

for lIdx=1:length(list)
    dataSfname = sprintf('neura-sparse01/imu/%s', list{lIdx}.file);
    load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, list{lIdx}.file));
    
    targetname = sprintf('%s/%s-%s-%s', targetdir, ns, list{lIdx}.file, list{lIdx}.name);
    
    %% algorithm 1
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo1));
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
    vel = vb.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
    sensors.PELVVelRef = vel.MIDPEL;
    sensors.LANKVelRef = vel.LTIO;
    sensors.RANKVelRef = vel.RTIO;
    sensors = experiment.buildSensorStructFromDebug(sensors, estState, estState2, cs.est);
    
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
    elseif strcmp(cs.stepDetection, 'av03')
        csGfrAcc = gfrAcc.(eLabel);
        
        csNSamples = size(csGfrAcc.MP, 1);
        bIsStatMP = false(csNSamples, 1);
        bIsStatLA = revStepDetect.stepL(idx);
        bIsStatRA = revStepDetect.stepR(idx);
    end
   
    % side by side
%     estBody1S2S = estBody1Rel.toWorldFrame(vb.MIDPEL, vb.qRPV);
    estBody1S2S = estBody1Rel.toWorldFrame(vb.MIDPEL, estBody1.qRPV);
    viconBodyS2S = viconBodyRel.toWorldFrame(vb.MIDPEL+[0.5 0 0], vb.qRPV);
    
    [estBody1Debug, sensorsDebug] = experiment.buildgrBodyDebug(estBody1, sensors, estState2, cs.est);
    [viconBodyDebug, ~] = experiment.buildgrBodyDebug(vb, struct(), false, 'vicon');
    bIsStatLADebug = repelem(bIsStatLA, 3, 1);
    bIsStatRADebug = repelem(bIsStatRA, 3, 1);
    
    % velocity constraint debug
    sensors = experiment.buildVelCstrDebug(sensors, estBody1, estState, cs.est);
    estState1Debug = experiment.buildStateDebug(estState, estState2, cs.est);
    sensorsDebug = experiment.buildVelCstrDebug(sensorsDebug, estBody1Debug, estState1Debug, cs.est);
    
    %% algorithm 2
    load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, list{lIdx}.file, list{lIdx}.algo2));
    
    estBody2 = estBody.copy();
    estBody2Rel = estBody2.changeRefFrame('MIDPEL');
    estBody2S2S = estBody2Rel.toWorldFrame(vb.MIDPEL+[0.5 0 0], estBody2.qRPV);
    [estBody2Debug, sensorsDebug2] = experiment.buildgrBodyDebug(estBody2, struct(), estState2, cs.est, 'Ref');
    
    sensors = experiment.buildSensorStructFromDebug(sensors, estState, estState2, cs.est, 'Ref');
    
    fn = fieldnames(sensorsDebug2);
    for i=1:length(fn)
        n = fn{i};
        sensorsDebug.(n) = sensorsDebug2.(n);
    end
    
    % velocity constraint debug
    sensors = experiment.buildVelCstrDebug(sensors, estBody2, estState, cs.est, 'Ref');
    estState2Debug = experiment.buildStateDebug(estState, estState2, cs.est);
    sensorsDebug = experiment.buildVelCstrDebug(sensorsDebug, estBody2Debug, estState2Debug, cs.est, 'Ref');
    
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