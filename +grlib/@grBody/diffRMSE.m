% ======================================================================
%> @brief Returns the RMSE difference between grBody1 and grBody2
%>
%> @param obj1 grBody 1 (self)
%> @param obj1 grBody 2 (other)
%> @param ref reference. default: MIDPEL
%> @param seq orientation sequence. default: YXZ
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = diffRelRMSE(obj1, obj2, ref, seq)
    if nargin <= 2
        ref = 'MIDPEL';
        seq = 'YXZ';
    end
    
    rawDiff = obj1.diff(obj2, seq);
    posFields = obj1.posList;
    
    out = struct;
    out.posMeanRMSE = 0.0;
    out.oriMeanRMSE = 0.0;
    
    val = [];
    for i=1:length(posFields)
        d = rawDiff.(posFields{i});
        val(end+1,:) = sqrt(nanmean(d.^2, 1));
        out.(posFields{i}) = val(end,:);
    end
    out.posMeanRMSE = mean(val(:));
    
    oriFields = obj1.oriList;
    val = [];
    for i=1:length(oriFields)
        d = rawDiff.(oriFields{i});
        val(end+1,:) = sqrt(nanmean(d.^2, 1));
        out.(oriFields{i}) = val(end,:);
    end
end