function obj = loadCalib(fname)
	% Load calibration files (calib_*.txt)
	%
	% :param fname: .sensors filename
	%
	% :return: output TCDBody with calibration data
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    %% Check function input
    validateattributes(fname, {'char'}, {});

    %% Variable initialization
    obj = mocapdb.XsensBody('srcFileName', fname, 'nSamples', 1, ...
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