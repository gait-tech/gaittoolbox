function out = loadCSV(fname, idxIfRaw)
	% Load CSV files
	%
	% :param fname: input filename
    % :param idxIfRaw: subset index if file contains untrimmed distance measurements
    % :type idxIfRaw: length 2 array of first and last index
	%
	% :return: out - table containing distances
	%
	% .. Author: - Luke Sy (UNSW GSBME) 17 Jun 2020
    
    out = readtable(fname);
    if nargin <= 1
        idxIfRaw = [1, size(out,1)];
    end
    fid = fopen(fname, 'r');
    buf = sscanf(fgetl(fid), '// IsRaw=%d Version=%f\n');
    fclose(fid);
    if idxIfRaw(2) > size(out, 1)
        idxIfRaw(2) = size(out, 1);
    end
    if isempty(buf) || buf(1)
        out = out(idxIfRaw(1):idxIfRaw(2), :);
    end
end

