% ======================================================================
%> @brief Get end index
%>
%> @retval endIdx last index
% ======================================================================
function endIdx = getEndIndex(obj)
    endIdx = length(obj.Hips(:,1));
end