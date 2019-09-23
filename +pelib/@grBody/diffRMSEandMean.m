% ======================================================================
%> @brief Returns the RMSE and Mean difference between grBody1 and grBody2
%>
%> @param obj1 grBody 1 (self)
%> @param obj1 grBody 2 (other)
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = diffRMSEandMean(obj1, obj2)
    seq = 'YXZ';
    
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
            out.(sprintf('%sRMSEnobias', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sStd', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sMean', oriFields{i})) = [nan nan nan];
            out.(sprintf('%sCorrCoef', oriFields{i})) = [nan nan nan];
        end
        out.oriMeanRMSE = nan;
        out.oriMeanMean = nan;
        out.dOri = nan;
        out.dOrinobias = nan;
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

        buf = {}; objList = {obj1, obj2};
        for i=1:2
            buf{i} = {};
            [r1,r2,r3] = quat2angle(objList{i}.qRPV, seq);
            buf{i}.qRPV  = [r2 r1 r3]*180/pi;
            buf{i}.qLHIP = objList{i}.calcJointAnglesLHip()*180/pi;
            buf{i}.qRHIP = objList{i}.calcJointAnglesRHip()*180/pi;
            buf{i}.qLKNE = objList{i}.calcJointAnglesLKnee()*180/pi;
            buf{i}.qRKNE = objList{i}.calcJointAnglesRKnee()*180/pi;
        end      
        nOriFields = length(oriFields);
        valRMSE = zeros(nOriFields, 3);
        valRMSEnobias = zeros(nOriFields, 3);
        valMean = zeros(nOriFields, 3);
        valStd =  zeros(nOriFields, 3);
        valCorrCoef = zeros(nOriFields, 3);
        for i=1:nOriFields
            d = rawDiff.(oriFields{i});
            valMean(i,:) = nanmean(d, 1);
            valRMSE(i,:) = sqrt(nanmean(d.^2, 1));
            valRMSEnobias(i,:) = sqrt(nanmean((d-valMean(i,:)).^2, 1));
            valStd(i,:) = nanstd(d, 1);
            theta1 = buf{1}.(oriFields{i});
            theta2 = buf{2}.(oriFields{i});
            thetaIdx = ~isnan(theta1) & ~isnan(theta2);
            R = zeros(2,2,3);
            for j=1:3
                R(:,:,j) = corrcoef(theta1(thetaIdx(:,j),j), ...
                                    theta2(thetaIdx(:,j),j));
            end
            valCorrCoef(i,:) = squeeze(R(2,1,:));
            out.(sprintf('%sRMSE', oriFields{i})) = valRMSE(i, :);
            out.(sprintf('%sRMSEnobias', oriFields{i})) = valRMSEnobias(i, :);
            out.(sprintf('%sStd',  oriFields{i})) = valStd(i, :);
            out.(sprintf('%sMean', oriFields{i})) = valMean(i,:);
            out.(sprintf('%sCorrCoef', oriFields{i})) = valCorrCoef(i,:);
        end
        out.oriMeanRMSE = mean(valRMSE(:));
        out.oriMeanMean = mean(valMean(:));
        
        out.dOri = nanmean(obj1.calcDOri(obj2));
        out.dOrinobias = nanmean(obj1.calcDOrinobias(obj2));
        out.dPos = nanmean(obj1.calcDPos(obj2));
    end
    out.trialDuration = obj1.nSamples;
end