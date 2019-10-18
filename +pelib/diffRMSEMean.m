% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
% @brief Returns the mean of diffRMSE results
%
%
% :param rs array of diffRMSE results
%
% :return: out returns the mean of diffRMSE results
% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
function out = diffRMSEMean(rs)
    rs2 = struct2table(rs);
    
    varNames = rs2.Properties.VariableNames;
    
    out = struct;
    for i=1:length(varNames)
        out.(varNames{i}) = mean(rs2.(varNames{i}), 1);
    end
end