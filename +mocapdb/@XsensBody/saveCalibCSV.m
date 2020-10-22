function saveCalibCSV(obj, fname, mode)
	% Save calibration as CSV file
	%
	% :param obj: this XsensBody
    % :type obj: :class:`+mocapdb.@XsensBody`
	% :param fname: filename of file to be saved
	% :param mode: save calib mode (bit 1 = qOri, 2 = acc_bias, 3 = both)
    % :type mode: Optional, integer. Defaults to 1.
    %
	% .. Author: - Luke Sy (UNSW GSBME)

    if nargin <= 2
        mode = 1;
    end
    
    t = table;
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if sum(size(obj.(n))) ~= 0
            if mode == 1
                t = [t; {n, obj.(n).ori(1,:)}];
            elseif mode == 2
                t = [t; {n, obj.(n).acc(1,:)}];
            else
                t = [t; {n, obj.(n).ori(1,:), obj.(n).acc(1,:)}];
            end
        end
    end
    if mode == 1
        t.Properties.VariableNames = {'name', 'q_ori'};
    elseif mode == 2
        t.Properties.VariableNames = {'name', 'acc'};
    else
        t.Properties.VariableNames = {'name', 'q_ori', 'acc'};
    end
    writetable(t, fname);
end