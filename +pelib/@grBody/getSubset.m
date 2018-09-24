% ======================================================================
%> @brief Get subset of grBody
%> @author Luke Sy (UNSW GSBME)
%> @date 22 Sept 2018
%>
%> Example:
%>      out = obj.getSubset(5:100, segAlias);
%>
%> @param obj class grBody (self)
%> @param idx indices of data to be included in out
%> @retval out grBody class whose data only includes the rows in idx
% ======================================================================
function out = getSubset(obj, idx)
    validateattributes(idx, {'numeric'}, {});
    
    out = obj.copy();
    segList = [obj.posList obj.oriList];
    
    for i=1:length(segList)
        n = segList{i};
        if sum(size(out.(n))) ~= 0
            out.(n) = out.(n)(idx, :);
        end
    end
    out.nSamples = length(idx);
end