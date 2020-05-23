function [obj, idx] = loadCSV(fname)
	% Load CSV files (opposite of exportCSV)
	%
	% :param fname: input file
	%
	% :return: obj - ViconBody
	%
	% .. Author: - Luke Sy (UNSW GSBME) 23 May 2020

    obj = mocapdb.ViconBody('srcFileName', fname);
    if ~endsWith(fname, '.csv')
        fname = sprintf("%s.csv", fname);
    end
    
    T = readtable(fname);
    for i = obj.posList
        cn = [sprintf("%s_1", i{1}), sprintf("%s_2", i{1}), ...
              sprintf("%s_3", i{1})];
        % if attribute is a column in table T
        if any(strcmp(cn(1), T.Properties.VariableNames))
            obj.(i{1}) = [T.(cn(1)) T.(cn(2)) T.(cn(3))];
        end
    end
    for i = obj.oriList
        cn = [sprintf("%s_1", i{1}), sprintf("%s_2", i{1}), ...
              sprintf("%s_3", i{1}), sprintf("%s_4", i{1})];
        % if attribute is a column in table T
        if any(strcmp(cn(1), T.Properties.VariableNames))
            obj.(i{1}) = [T.(cn(1)) T.(cn(2)) T.(cn(3)) T.(cn(4))];
        end
    end
    fid = fopen(fname, 'r');
    buf = sscanf(fgetl(fid), '// fs=%d nSamples=%d\n');
    obj.fs = buf(1); obj.nSamples = buf(2);
    obj.frame = sscanf(fgetl(fid), '// frame=%s\n');
    obj.posUnit = sscanf(fgetl(fid), '// posUnit=%s\n');
    fgetl(fid);
    idx = sscanf(fgetl(fid), '// idx=%d,%d\n')';
    fclose(fid);
end