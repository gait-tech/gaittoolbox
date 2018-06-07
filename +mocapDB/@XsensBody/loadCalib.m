% ======================================================================
%> @brief Load calib_*.txt file
%>
%> Load the calibration file
%>
%> @param fname .sensors filename
%>
%> @retval output TCDBody with calibration data
% ======================================================================
function obj = loadCalib(fname)
    %% Check function input
    validateattributes(fname, {'char'}, {});

    %% Variable initialization
    obj = tcdlib.XsensBody('srcFileName', fname, 'nSamples', 1, ...
                           'frame', 'calib');
    
    %% Load calibration data
    fileID = fopen(fname, 'r');
    colN = fscanf(fileID, '%d', 1);
    
    for i=1:colN
        rname = fscanf(fileID, '%s', 1);
        rval = fscanf(fileID, '%f', [1,4]);
        rval = [rval(4) rval(1:3)];
        obj.(rname) = table(rval, 'VariableNames', {'ori'});
    end
    fclose(fileID);
end