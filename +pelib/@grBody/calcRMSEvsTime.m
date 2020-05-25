function out = calcRMSEvsTime(obj1, obj2, includeRoot, targetSeg)
	% Returns the RMSE difference between grBody1 and grBody2 VS time
	%
	% :param obj1: grBody 1 (self)
	% :param obj1: grBody 2 (other)
    % :param includeRoot: [boolean] if root/pelvis is included
    % :param targetSeg: segments to be computed (usually occluded)
	%
	% :return: out - struct with the difference of pos and ori parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME) 2020/04/09
    
    if nargin <= 2, includeRoot = true; end
    if nargin <= 3, targetSeg = {'qLTH', 'qRTH'}; end
    seq = 'YXZ';
    
    out = struct;
    posFields = obj1.posList;
    oriFields = {'qRPV', 'qLHIP', 'qRHIP', 'qLKNE', 'qRKNE', 'qLANK', 'qRANK'};
    
    if ~isobject(obj2) && isnan(obj2)
        for i=1:length(posFields)
            out.(sprintf('%sAbsErr', posFields{i})) = [nan nan nan];
            out.(sprintf('%sErr', posFields{i})) = [nan nan nan];
        end
        for i=1:length(oriFields)
            out.(sprintf('%sAbsErr', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sAbsErrnobias', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sErr', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sErrnobias', oriFields{i})) = [nan nan nan];
        end
        out.dOri = nan;
        out.dOrinobias = nan;
        out.dPos = nan;
    else
        rawDiff = obj1.diff(obj2, seq);
        
        nPosFields = length(posFields);
        for i=1:nPosFields
            d = rawDiff.(posFields{i});
            if isempty(d)
                d = nan(size(rawDiff.MIDPEL));
            end
            
            out.(sprintf('%sAbsErr', posFields{i})) = abs(d);
            out.(sprintf('%sErr', posFields{i})) = d;
        end
        % summary for joint positions
        
        buf = {}; objList = {obj1, obj2};
        for i=1:2
            buf{i} = {};
            [r1,r2,r3] = quat2angle(objList{i}.qRPV, seq);
            buf{i}.qRPV  = [r2 r1 r3]*180/pi;
            buf{i}.qLHIP = objList{i}.calcJointAnglesLHip()*180/pi;
            buf{i}.qRHIP = objList{i}.calcJointAnglesRHip()*180/pi;
            buf{i}.qLKNE = objList{i}.calcJointAnglesLKnee()*180/pi;
            buf{i}.qRKNE = objList{i}.calcJointAnglesRKnee()*180/pi;
            buf{i}.qLANK = objList{i}.calcJointAnglesLAnkle()*180/pi;
            buf{i}.qRANK = objList{i}.calcJointAnglesRAnkle()*180/pi;
        end      
        nOriFields = length(oriFields);
        valMean = struct();
        for i=1:nOriFields
            d = rawDiff.(oriFields{i});
            if isempty(d)
                bufMean = [nan nan nan];
                d = nan(size(rawDiff.MIDPEL));
            else
                bufMean = nanmean(d, 1);
            end
            valMean.(oriFields{i}) = bufMean;
            out.(sprintf('%sAbsErr', oriFields{i})) = abs(d);
            out.(sprintf('%sAbsErrnobias', oriFields{i})) = abs(d-bufMean);
            out.(sprintf('%sErr', oriFields{i})) = d;
            out.(sprintf('%sErrnobias', oriFields{i})) = d-bufMean;
        end
        
        % summary for joint angles
        d = [rawDiff.qLHIP rawDiff.qRHIP rawDiff.qLKNE(:,2) rawDiff.qRKNE(:,2)];
        out.qHipKneeYRMSE = sqrt(nanmean(d.^2,2));
        d = cat(3,rawDiff.qLHIP,rawDiff.qRHIP);
        out.qHipRMSE = sqrt(nanmean(d.^2,3));
        d = cat(3,rawDiff.qLKNE,rawDiff.qRKNE);
        out.qKneeRMSE = sqrt(nanmean(d.^2,3));
        d = cat(3,rawDiff.qLANK,rawDiff.qRANK);
        out.qAnkleRMSE = sqrt(nanmean(d.^2,3));
        
        d = [rawDiff.qLHIP-valMean.qLHIP rawDiff.qRHIP-valMean.qRHIP ...
             rawDiff.qLKNE(:,2)-valMean.qLKNE(2) ...
             rawDiff.qRKNE(:,2)-valMean.qRKNE(2)];
        out.qHipKneeYRMSEnobias = sqrt(nanmean(d.^2,2));
        d = cat(3,rawDiff.qLHIP-valMean.qLHIP,rawDiff.qRHIP-valMean.qRHIP);
        out.qHipRMSEnobias = sqrt(nanmean(d.^2,3));
        d = cat(3,rawDiff.qLKNE-valMean.qLKNE,rawDiff.qRKNE-valMean.qRKNE);
        out.qKneeRMSEnobias = sqrt(nanmean(d.^2,3));
        d = cat(3,rawDiff.qLANK-valMean.qLANK,rawDiff.qRANK-valMean.qRANK);
        out.qAnkleRMSEnobias = sqrt(nanmean(d.^2,3));
        
        out.dOri = obj1.calcDOri(obj2,targetSeg);
        [out.dOrinobias, out.dOribias] = obj1.calcDOrinobias(obj2,targetSeg);
        out.dOribias = repmat(out.dOribias, size(out.dOri));
        out.dPos = obj1.calcDPos(obj2,includeRoot);
    end
%     out.trialDuration = obj1.nSamples;
end