function obj = loadCalibCSV(fname)
	% Save calibration as CSV file
	%
	% :param obj: this XsensBody
	% :param fname: filename of file to be saved
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    T = readtable(fname);
    
    obj = mocapdb.XsensBody();
    for i=1:size(T, 1)
        n = T.name{i};
        obj.(n) = table;
        if ismember('q_ori_1', T.Properties.VariableNames)
            obj.(n).ori = T{i,{'q_ori_1', 'q_ori_2', 'q_ori_3', 'q_ori_4'}};
        end
        if ismember('acc_1', T.Properties.VariableNames)
            obj.(n).acc = T{i,{'acc_1', 'acc_2', 'acc_3'}};
        end
    end
end