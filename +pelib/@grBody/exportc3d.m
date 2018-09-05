% ======================================================================
%> @brief generate c3d file
%>
%> @param obj grBody (self)
%> @param fname output file name
%>
%> @retval acq pointer to new btk c3d
% ======================================================================
function acq = exportc3d(obj, fname)
    acq = btkNewAcquisition(length(obj.posList), obj.nSamples);
    btkSetFrequency(acq, obj.fs);
    
    % set points
    for i = 1:length(obj.posList)
        ptName = obj.posList{i};
        btkSetPointLabel(acq, i, ptName)
        btkSetPoint(acq, i, (obj.(ptName))*1000) % values: array of 20 by 3 by 2000 
    end
    btkWriteAcquisition(acq, fname);
end