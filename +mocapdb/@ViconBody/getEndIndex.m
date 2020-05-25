<<<<<<< HEAD
function endIdx = getEndIndex(obj)
	% Get end index
	%
	% :return: endIdx - last index
	%
	% .. Author: - Luke Sy (UNSW GSBME)
    endIdx = length(obj.PELV(:,1));
=======
function endIdx = getEndIndex(obj, untilfirstnan)
	% Get end index
	%
    % :param obj: this ViconBody
    % :param untilfirstnan: default False. if True, return end index at
    % first Nan
    %
	% :return: endIdx - last index
	%
	% .. Author: - Luke Sy (UNSW GSBME)    
    if nargin <= 1
        untilfirstnan = false;
    end
    
    endIdx = length(obj.PELV(:,1));
    if untilfirstnan
        startIdx = obj.getStartIndex();
        for i=1:length(obj.posList)
            buf = obj.(obj.posList{i});
            if ~isempty(buf)
                idx = find(any(isnan(buf(startIdx:end,:)), 2), 1);
                if ~isempty(idx)
                    endIdx = min(idx(1)-2+startIdx, endIdx);
                end
            end
        end
        for i=1:length(obj.oriList)
            buf = obj.(obj.oriList{i});
            if ~isempty(buf)
                idx = find(any(isnan(buf(startIdx:end,:)), 2), 1);
                if ~isempty(idx)
                    endIdx = min(idx(1)-2+startIdx, endIdx);
                end
            end
        end
    end
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
end