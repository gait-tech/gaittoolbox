% ======================================================================
%> Generate step detect templates from variance
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
stepDir = sprintf('%s/step-detect', dir);
outDir = sprintf('%s/step-detect-c3d', dir);
ns = "NS2";
algo = "NS2+Aw__sOw__sIw__v+Sav01+M76+C355";
outEdit = sprintf('%s/edit-v3.csv', dir);

velThreshold = 0.1;
minStepFrameLength = 5;

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

    name = sprintf("%s-%s", n.subj, n.act);
    imuStepFName = sprintf("%s-imuStepDetect.csv", name);
    revStepFName = sprintf("%s-revStepDetect.csv", name);
    imuStep = readtable(sprintf("%s/%s", stepDir, imuStepFName));
    viconStep = readtable(sprintf('%s/%s-viconStepDetect.csv', stepDir, name));
    load(sprintf('%s/%s-%s-debug.mat', expDir, ns, name));
    load(sprintf('%s/%s-%s-%s.mat', expDir, ns, name, algo));
    
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
    
    bIsStatLA = imuStep.stepL(idx);
    bIsStatRA = imuStep.stepR(idx);    
%     bIsStatLA = viconStep.stepL(idx);
%     bIsStatRA = viconStep.stepR(idx);    
    nSamples = length(idx);
    speventsL = false(nSamples, 1);
    speventsR = false(nSamples, 1);
    
    outInst(outIdx) = struct('cmd', 'load', 'v0', imuStepFName, 'v1', "", 'v2', "", 'cmt',  "");
    outIdx = outIdx+1;
    %% check those that are detected as double step from IMU data but has some velocity according from Vicon data
    velmagLA = vecnorm(vel.LTIO, 2, 2);    
    velmagRA = vecnorm(vel.RTIO, 2, 2);
    doublefoot = (bIsStatLA & bIsStatRA) & ((velmagLA > velThreshold) | (velmagRA > velThreshold));
    [sIdxs eIdxs] = getStartEndIndices([0; doublefoot; 0]);
    for j=1:length(sIdxs)
        sIdx = sIdxs(j)-1; eIdx = min(eIdxs(j)-1, nSamples);
        velmag2LA = norm(velmagLA(sIdx:eIdx));
        velmag2RA = norm(velmagRA(sIdx:eIdx));
        if velmag2LA > velmag2RA, side = "L";
        elseif velmag2LA < velmag2RA, side = "R";
        else side = ""; end

        outInst(outIdx) = struct('cmd', 'clear', 'v0', side, 'v1', sIdx, 'v2', eIdx, 'cmt',  "doublefoot");
        outIdx = outIdx+1;
    end
    speventsL = speventsL | doublefoot;
    speventsR = speventsR | doublefoot;

    velLA = ~speventsL & bIsStatLA & (velmagLA > velThreshold);
    [sIdxs eIdxs] = getStartEndIndices([0; velLA; 0]);
    for j=1:length(sIdxs)
        sIdx = sIdxs(j)-1; eIdx = min(eIdxs(j)-1, nSamples);
        outInst(outIdx) = struct('cmd', 'clear', 'v0', "L", 'v1', sIdx, 'v2', eIdx, 'cmt',  "velThreshL");
        outIdx = outIdx+1;
    end
    speventsL = speventsL | velLA;

    shortLA = ~speventsL & bIsStatLA;
    [sIdxs eIdxs] = getStartEndIndices([0; shortLA; 0]);
    for j=1:length(sIdxs)
        sIdx = sIdxs(j)-1; eIdx = min(eIdxs(j)-1, nSamples);
        if ((eIdx - sIdx) <= minStepFrameLength) && (eIdx ~= nSamples)
            outInst(outIdx) = struct('cmd', 'clear', 'v0', "L", 'v1', sIdx, 'v2', eIdx, 'cmt',  "shortFrameThreshL");
            outIdx = outIdx+1;
        else
            shortLA(sIdx:eIdx) = false;
        end
    end
    speventsL = speventsL | shortLA;

    velRA = ~speventsR & bIsStatRA & (velmagRA > velThreshold);
    [sIdxs eIdxs] = getStartEndIndices([0; velRA; 0]);
    for j=1:length(sIdxs)
        sIdx = sIdxs(j)-1; eIdx = min(eIdxs(j)-1, nSamples);
        outInst(outIdx) = struct('cmd', 'clear', 'v0', "R", 'v1', sIdx, 'v2', eIdx, 'cmt',  "velThreshR");
        outIdx = outIdx+1;
    end
    speventsR = speventsR | velRA;

    shortRA = ~speventsR & bIsStatRA;
    [sIdxs eIdxs] = getStartEndIndices([0; shortRA; 0]);
    for j=1:length(sIdxs)
        sIdx = sIdxs(j)-1; eIdx = min(eIdxs(j)-1, nSamples);
        if ((eIdx - sIdx) <= minStepFrameLength) && (eIdx ~= nSamples)
            outInst(outIdx) = struct('cmd', 'clear', 'v0', "R", 'v1', sIdx, 'v2', eIdx, 'cmt', "shortFrameThreshR");
            outIdx = outIdx+1;
        else
            shortRA(sIdx:eIdx) = false;
        end
    end
    speventsR = speventsR | shortRA;

    outInst(outIdx) = struct('cmd', 'save', 'v0', revStepFName, 'v1', '', 'v2', '', 'cmt', "");
    outIdx = outIdx+1;
    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, struct(), 1, speventsL|speventsR);
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
end
writetable(struct2table(outInst), outEdit);

function [sIdx eIdx] = getStartEndIndices(events)
    events2 = events - events([1, 1:end-1]);
    [sIdx, val] = find(events2 == 1);
    [eIdx, val] = find(events2 == -1);
end