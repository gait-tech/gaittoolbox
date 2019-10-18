function startIdx = getStartIndex(obj)
	% Get first valid index
	%
	% :return: startIdx - first valid index
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    startIdx = 1;
    for i=1:length(obj.posList)
        idx = find(~any(isnan(obj.(obj.posList{i})), 2), 1);
        startIdx = max(idx(1), startIdx);
    end
    for i=1:length(obj.oriList)
        idx = find(~any(isnan(obj.(obj.oriList{i})), 2), 1);
        startIdx = max(idx(1), startIdx);
    end
end