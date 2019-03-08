% ======================================================================
%> Generate step detect templates from variance
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
stepDir = sprintf('%s/step-detect', dir);
% outDir = sprintf('%s/step-detect-c3d', dir);
outDir = 'C:\Users\z5151460\OneDrive - UNSW\Thesis - Sparse Mocap\Aim 1\Step Detection Fix\step-detect-v3-c3d';
ns = "NS2";
algo = "NS2+Aw__sOw__sIw__v+Sav01+M76+C355";
outEdit = sprintf('%s/edit-v3dynamic.csv', dir);

stepFrameLength = 5;
stepHalfFrameLength = idivide(int32(stepFrameLength), 2, 'floor');

dataList = readtable(sprintf('%s/data-list-v2.csv', dir));
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataN = size(dataList, 1);

clear outInst;
outIdx = 1;

for i = 1:dataN
    n = table2struct(dataList(i, :));
    if isempty(strfind(n.act, "HighKneeJog")) && ...
       isempty(strfind(n.act, "Jog")) && ...
       isempty(strfind(n.act, "JumpingJacks"))
        continue;
    end
    
    name = sprintf("%s-%s", n.subj, n.act);
    imuStepFName = sprintf("%s-imuStepDetect.csv", name);
    revStepFName = sprintf("%s-revStepDetect.csv", name);
    imuStep = readtable(sprintf("%s/%s", stepDir, imuStepFName));
    viconStep = readtable(sprintf('%s/%s-viconStepDetect.csv', stepDir, name));
    load(sprintf('%s/%s-%s-debug.mat', expDir, ns, name));
    load(sprintf('%s/%s-%s-%s.mat', expDir, ns, name, algo), 'cs', 'estBody', 'estState');
    
    targetname = sprintf('%s/%s-%s-%s', outDir, ns, name, algo);
    
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
    nSamples = length(idx);
    speventsL = false(nSamples, 1);
    speventsR = false(nSamples, 1);
    
    sensors = struct();
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
    
    if strcmp(cs.stepDetection, 'av01')
        bIsStatLA = imuStep.stepL(idx);
        bIsStatRA = imuStep.stepR(idx);  
    elseif strcmp(cs.stepDetection, 'av02')
        bIsStatLA = viconStep.stepL(idx);
        bIsStatRA = viconStep.stepR(idx);
    else
        display('step detection style not supported');
        bIsStatLA = false(nSamples, 1);
        bIsStatRA = false(nSamples, 1);
    end
    
    outInst(outIdx) = struct('cmd', 'load', 'v0', imuStepFName, 'v1', "", 'v2', "", 'cmt',  "");
    outIdx = outIdx+1;
    
    %% find all lowest points
    [pks, locs] = findpeaks(-vb.LTIO(:,3));
    for j=1:length(locs)
        sIdx = max(1, locs(j)-stepHalfFrameLength);
        eIdx = min(locs(j)+stepHalfFrameLength, nSamples);
        speventsL(sIdx:eIdx) = true;
        
        outInst(outIdx) = struct('cmd', 'set', 'v0', "L", 'v1', sIdx, 'v2', eIdx, 'cmt',  "lowestpoint");
        outIdx = outIdx+1;
    end
    
    [pks, locs] = findpeaks(-vb.RTIO(:,3));
    for j=1:length(locs)
        sIdx = max(1, locs(j)-stepHalfFrameLength);
        eIdx = min(locs(j)+stepHalfFrameLength, nSamples);
        speventsR(sIdx:eIdx) = true;
        
        outInst(outIdx) = struct('cmd', 'set', 'v0', "R", 'v1', sIdx, 'v2', eIdx, 'cmt',  "lowestpoint");
        outIdx = outIdx+1;
    end

    outInst(outIdx) = struct('cmd', 'save', 'v0', revStepFName, 'v1', '', 'v2', '', 'cmt', "");
    outIdx = outIdx+1;
    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, struct(), 1, speventsL|speventsR);
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
end
writetable(struct2table(outInst), outEdit);