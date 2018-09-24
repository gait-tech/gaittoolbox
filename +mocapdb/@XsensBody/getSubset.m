% ======================================================================
%> @brief Get subset of xsens measurements
%> @author Luke Sy (UNSW GSBME)
%> @date 20 Sept 2018
%>
%> Example:
%>      out = obj.getSubset(5:100, segAlias);
%>
%> @param obj class XsensBody (self)
%> @param idx indices of data to be included in out
%> @retval out XsensBody class whose data only includes the rows in idx
% ======================================================================
function out = getSubset(obj, idx)
    validateattributes(idx, {'numeric'}, {});
    
    out = obj.copy();
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if sum(size(out.(n))) ~= 0
            out.(n) = out.(n)(idx, :);
        end
    end
    out.nSamples = length(idx);
end