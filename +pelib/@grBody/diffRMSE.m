function out = diffRMSE(obj1, obj2, ref, seq)
	% Returns the RMSE difference between grBody1 and grBody2
	%
	% :param obj1: grBody 1 (self)
	% :param obj1: grBody 2 (other)
	% :param ref: reference. default: MIDPEL
	% :param seq: orientation sequence. default: YXZ
	%
	% :return: out - struct with the difference of pos and ori parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    if nargin <= 2
        ref = 'MIDPEL';
        seq = 'YXZ';
    end
    
    out = struct;
    posFields = obj1.posList;
%     oriFields = obj1.oriList;
    oriFields = {'qRPV', 'qLHIP', 'qRHIP', 'qLKNE', 'qRKNE'};
    
    if ~isobject(obj2) && isnan(obj2)
        for i=1:length(posFields)
            out.(posFields{i}) = [nan nan nan];
            out.(sprintf('%sStd', posFields{i})) = [nan nan nan];
        end
        out.posMeanRMSE = nan;
        
        for i=1:length(oriFields)
            out.(oriFields{i}) = [nan nan nan];
            out.(sprintf('%sStd', oriFields{i})) = [nan nan nan];
        end
        out.oriMeanRMSE = nan;
        out.dOri = nan;
        out.dPos = nan;
    else
        rawDiff = obj1.diff(obj2, seq);

        out.posMeanRMSE = 0.0;
        out.oriMeanRMSE = 0.0;

        val = [];
        valStd = [];
        for i=1:length(posFields)
            d = rawDiff.(posFields{i});
            if(~isempty(d))
                val(end+1,:) = sqrt(nanmean(d.^2, 1));
                valStd(end+1,:) = nanstd(d, 1);
            else
                val(end+1,:) = [nan nan nan];
                valStd(end+1,:) = [nan nan nan];
            end
            out.(posFields{i}) = val(end,:);
            out.(sprintf('%sStd', posFields{i})) = valStd(end,:);
        end
        out.posMeanRMSE = nanmean(val(:));

        val = [];
        valStd = [];
        for i=1:length(oriFields)
            d = rawDiff.(oriFields{i});
            if(~isempty(d))
                val(end+1,:) = sqrt(nanmean(d.^2, 1));
                valStd(end+1,:) = nanstd(d, 1);
            else
                val(end+1,:) = [nan nan nan];
                valStd(end+1,:) = [nan nan nan];
            end
            out.(oriFields{i}) = val(end,:);
            out.(sprintf('%sStd', oriFields{i})) = valStd(end,:);
        end
        out.oriMeanRMSE = nanmean(val(:));
        
        out.dOri = nanmean(obj1.calcDOri(obj2));
        out.dPos = nanmean(obj1.calcDPos(obj2));
    end
    out.trialDuration = obj1.nSamples;
end