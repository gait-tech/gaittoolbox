function out = getMean(obj)
	% Calculate the mean
	% 
	% :param obj: this object
    % :type obj: :class:`+mocapdb\@XsensBody`
	%
	% :return: out - XsensBody with 1 row per attribute containing mean data
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 Jun 9
    
    %% Variable initialization
    out = mocapdb.XsensBody('srcFileName', obj.srcFileName, ...
                            'nSamples', 1, ...
                            'frame', 'calib');
    
    %% Calculation
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if isempty(obj.(n)), continue; end;
        for j=["acc", "gyr", "mag"]
            if ~isempty(obj.(n).(j))
                out.(n).(j) = mean(obj.(n).(j), 1);
            end
        end
        q = quaternion(obj.(n).ori);
        if ~isempty(q)
            out.(n).ori = compact(meanrot(q, 'omitnan'));
        end
    end
end