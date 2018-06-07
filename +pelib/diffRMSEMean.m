% ======================================================================
%> @brief Returns the mean of diffRMSE results
%>
%>
%> @param rs array of diffRMSE results
%>
%> @retval out returns the mean of diffRMSE results
% ======================================================================
function out = diffRMSEMean(rs)
    rs2 = struct2table(rs);
    
    varNames = rs2.Properties.VariableNames;
    
    out = struct;
    for i=1:length(varNames)
        out.(varNames{i}) = mean(rs2.(varNames{i}), 1);
    end
end