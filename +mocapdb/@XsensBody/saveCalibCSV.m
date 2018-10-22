% ======================================================================
%> @brief Save calibration as CSV file
%>
%> @param obj this XsensBody
%> @param fname filename of file to be saved
% ======================================================================
function saveCalibCSV(obj, fname)
    t = table;
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if sum(size(obj.(n))) ~= 0
            t = [t; {n, obj.(n).ori(1,:)}];
        end
    end
    
    t.Properties.VariableNames = {'name', 'q_ori'};
    writetable(t, fname);
end