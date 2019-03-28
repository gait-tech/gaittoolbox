% ======================================================================
%> Generate step detect templates from variance
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
stepDir = sprintf('%s/step-detect', dir);
% outDir = sprintf('explore_output', dir);
outDir = 'C:\Users\z5151460\OneDrive - UNSW\Thesis - Sparse Mocap\Aim 1\Analysis - 20190312 - ckf results';
ns = "NS2";
% algo = "NS2+pfv1+Aw__sOw__sIw__v+Sav03+P001+M001";
algo = "";
velThreshold = 0.1;

dataList = readtable(sprintf('%s/data-list-v2.csv', dir));
options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
dataN = size(dataList, 1);

for algoB = {"NS2+Aw__vOw__vIw__v+Sav03+M00+C000"}% "NS2+Aw__sOw__sIw__v+Sav03+M76+C355"}
    algo = algoB{1};
    for i = 1:dataN
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
        idx = idx(1:estBody.nSamples);
        estBody = estBody.getSubset(1:estBody.nSamples);
        vb = vb.getSubset(1:estBody.nSamples);
        
        sensors = struct();
        sensors.PELVFreeAcc = gfrAcc.(eLabel).MP(1:estBody.nSamples,:);
        sensors.LANKFreeAcc = gfrAcc.(eLabel).LA(1:estBody.nSamples,:);
        sensors.RANKFreeAcc = gfrAcc.(eLabel).RA(1:estBody.nSamples,:);
        sensors.PELVFreeAccRef = gfrAcc.(aLabel).MP(1:estBody.nSamples,:);
        sensors.LANKFreeAccRef = gfrAcc.(aLabel).LA(1:estBody.nSamples,:);
        sensors.RANKFreeAccRef = gfrAcc.(aLabel).RA(1:estBody.nSamples,:);
        sensors.PELVVel = estState(1:estBody.nSamples, 4:6);
        sensors.LANKVel = estState(1:estBody.nSamples, 14:16);
        sensors.RANKVel = estState(1:estBody.nSamples, 24:26);
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
end

function [sIdx eIdx] = getStartEndIndices(events)
    events2 = events - events([1, 1:end-1]);
    [sIdx, val] = find(events2 == 1);
    [eIdx, val] = find(events2 == -1);
end