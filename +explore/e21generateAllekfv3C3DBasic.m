% ======================================================================
%> Generate step detect templates from variance
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
stepDir = sprintf('%s/step-detect', dir);
outDir = sprintf('explore_output', dir);
% outDir = 'C:\Users\z5151460\OneDrive - UNSW\Thesis - Sparse Mocap\Aim 1\Analysis - Step Fix - 01 Mar 2019\step-detect-v3-fixed-c3d';
ns = "NS2";
% algo = "NS2+pfv1+Aw__sOw__sIw__v+Sav03+P001+M001";
algo = "NS2+Aw__sOw__sIw__v+Sav03+M302+C000";
velThreshold = 0.1;

dataList = readtable(sprintf('%s/data-list-v2.csv', dir));
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataN = size(dataList, 1);

for i = [1 5 7 15]
    n = table2struct(dataList(i, :));
    
    name = sprintf("%s-%s", n.subj, n.act);
    imuStepFName = sprintf("%s-imuStepDetect.csv", name);
    revStepFName = sprintf("%s-revStepDetect.csv", name);
    imuStep = readtable(sprintf("%s/%s", stepDir, imuStepFName));
    viconStep = readtable(sprintf('%s/%s-viconStepDetect.csv', stepDir, name));
    load(sprintf('%s/%s-%s-debug.mat', expDir, ns, name));
    load(sprintf('%s/%s-%s-%s.mat', expDir, ns, name, algo));
    revStepDetect = readtable(sprintf('%s/step-detect/%s-%s-revStepDetect.csv', ...
                        dir, n.subj, n.act));

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
    uwbdist = estBody.calcMPLARAdist();
    sensors.MPLARAdist = [uwbdist.MPLA uwbdist.MPRA uwbdist.LARA];
    uwbdistRef = vb.calcMPLARAdist();
    sensors.MPLARAdistRef = [uwbdistRef.MPLA uwbdistRef.MPRA uwbdistRef.LARA];
    
    if strcmp(cs.stepDetection, 'av01')
        bIsStatLA = imuStep.stepL(idx);
        bIsStatRA = imuStep.stepR(idx);  
    elseif strcmp(cs.stepDetection, 'av02')
        bIsStatLA = viconStep.stepL(idx);
        bIsStatRA = viconStep.stepR(idx);    
    else % strcmp(cs.stepDetection, 'av03')
        bIsStatLA = revStepDetect.stepL(idx);
        bIsStatRA = revStepDetect.stepR(idx);
    end

    nSamples = length(idx);
    events = false(nSamples, 1);
    
    estBody.exportc3d(sprintf('%s.c3d', targetname), sensors, ...
                      vb, bIsStatLA, bIsStatRA, struct(), 1);
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
end

function [sIdx eIdx] = getStartEndIndices(events)
    events2 = events - events([1, 1:end-1]);
    [sIdx, val] = find(events2 == 1);
    [eIdx, val] = find(events2 == -1);
end