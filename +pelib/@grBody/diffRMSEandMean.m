% ======================================================================
%> @brief Returns the RMSE and Mean difference between grBody1 and grBody2
%>
%> @param obj1 grBody 1 (self)
%> @param obj1 grBody 2 (other)
%> @param seq orientation sequence. default: YXZ
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = diffRMSEandMean(obj1, obj2, seq)
    if nargin <= 2
        seq = 'YXZ';
    end
    
    out = struct;
    posFields = obj1.posList;
%     oriFields = obj1.oriList;
    oriFields = {'qRPV', 'qLHIP', 'qRHIP', 'qLKNE', 'qRKNE'};
    
    if ~isobject(obj2) && isnan(obj2)
        for i=1:length(posFields)
            out.(sprintf('%sRMSE', posFields{i})) = [nan nan nan];
            out.(sprintf('%sStd', posFields{i})) = [nan nan nan];
            out.(sprintf('%sMean', posFields{i})) = [nan nan nan];
        end
        out.posMeanRMSE = nan;
        out.posMeanMean = nan;
        
        for i=1:length(oriFields)
            out.(sprintf('%sRMSE', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sStd', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sMean', oriFields{i})) = [nan nan nan];
        end
        out.oriMeanRMSE = nan;
        out.oriMeanMean = nan;
        out.dOri = nan;
        out.dPos = nan;
    else
        rawDiff = obj1.diff(obj2, seq);

        out.posMeanRMSE = 0.0;
        out.posMeanMean = 0.0;
        out.oriMeanRMSE = 0.0;
        out.oriMeanMean = 0.0;
        
        nPosFields = length(posFields);
        valRMSE = zeros(nPosFields, 3);
        valMean = zeros(nPosFields, 3);
        valStd =  zeros(nPosFields, 3);
        for i=1:nPosFields
            d = rawDiff.(posFields{i});
            valMean(i,:) = nanmean(d, 1);
            valRMSE(i,:) = sqrt(nanmean(d.^2, 1));
            valStd(i,:) = nanstd(d, 1);
            out.(sprintf('%sRMSE', posFields{i})) = valRMSE(i, :);
            out.(sprintf('%sStd',  posFields{i})) = valStd(i, :);
            out.(sprintf('%sMean', posFields{i})) = valMean(i,:);
        end
        out.posMeanRMSE = mean(valRMSE(:));
        out.posMeanMean = mean(valMean(:));

        nOriFields = length(oriFields);
        valRMSE = zeros(nOriFields, 3);
        valMean = zeros(nOriFields, 3);
        valStd =  zeros(nOriFields, 3);
        for i=1:nOriFields
            d = rawDiff.(oriFields{i});
            valMean(i,:) = nanmean(d, 1);
            valRMSE(i,:) = sqrt(nanmean(d.^2, 1));
            valStd(i,:) = nanstd(d, 1);
            out.(sprintf('%sRMSE', oriFields{i})) = valRMSE(i, :);
            out.(sprintf('%sStd',  oriFields{i})) = valStd(i, :);
            out.(sprintf('%sMean', oriFields{i})) = valMean(i,:);
        end
        out.oriMeanRMSE = mean(valRMSE(:));
        out.oriMeanMean = mean(valMean(:));
        
        out.dOri = nanmean(obj1.calcDOri(obj2));
        out.dPos = nanmean(obj1.calcDPos(obj2));
    end
    out.trialDuration = obj1.nSamples;
end