% ======================================================================
%> @brief Returns the RMSE difference between grBody1 and grBody2
%>
%> @param obj1 grBody 1 (self)
%> @param obj1 grBody 2 (other)
%> @param seq orientation sequence
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = diffRMSE(obj1, obj2, seq)
    if nargin <= 2
        seq = 'YXZ';
    end
    
    rawDiff = obj1.diff(obj2, seq);
    fields = fieldnames(rawDiff);
    out = struct;
    
    for i=1:length(fields)
        d = rawDiff.(fields{i});
        out.(fields{i}) = sqrt(nanmean(d.^2, 1));
    end
end