% ======================================================================
%> @brief Get end index
%>
%> @retval endIdx last index
% ======================================================================
function endIdx = getEndIndex(obj)
    endIdx = length(obj.PELV(:,1));
end