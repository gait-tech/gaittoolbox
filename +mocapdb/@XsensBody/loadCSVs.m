function [obj, idx] = loadCSVs(fname)
	% Load CSV files (opposite of exportCSVs)
    % Pelvis, L_UpLeg, R_UpLeg, L_LowLeg, L_LowLeg2, 
    % R_LowLeg, R_LowLeg2, L_Foot, R_Foot
	% 
    % This function will specifically load <fname>-Pelvis.csv to obj.Pelvis,
    %   <fname>-L_Ankles to the obj.L_LowLeg, etc.
	%
	% :param fname: input filename
	%
	% :return: obj - XsensBody
	%
	% .. Author: - Luke Sy (UNSW GSBME) 23 May 2020

    % list of body segment names
    bs = mocapdb.XsensBody.segList;
    bsN = length(bs);
    
    obj = mocapdb.XsensBody('srcFileName', fname);
    for i = 1:bsN
        fpath = sprintf("%s-%s.csv", fname, bs{i});
        if exist(fpath)
            T = readtable(fpath);
            obj.(bs{i}) = table([T.ori_1 T.ori_2 T.ori_3 T.ori_4], ...
                                [T.acc_1 T.acc_2 T.acc_3], ...
                                [T.gyr_1 T.gyr_2 T.gyr_3], ...
                                [T.mag_1 T.mag_2 T.mag_3], ...
                               'VariableNames', {'ori', 'acc', 'gyr', 'mag'});
                           
            fid = fopen(fpath, 'r');
            buf = sscanf(fgetl(fid), '// fs=%d nSamples=%d\n');
            obj.fs = buf(1); obj.nSamples = buf(2);
            obj.frame = sscanf(fgetl(fid), '// frame=%s\n');
            fgetl(fid);
            idx = sscanf(fgetl(fid), '// idx=%d,%d\n')';
            fclose(fid);
        end
    end
end