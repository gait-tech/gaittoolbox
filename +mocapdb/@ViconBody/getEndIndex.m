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
        for i=1:length(obj.posList)
            idx = find(any(isnan(obj.(obj.posList{i})), 2), 1);
            if ~isempty(idx)
                endIdx = min(idx(1)-1, endIdx);
            end
        end
        for i=1:length(obj.oriList)
            idx = find(any(isnan(obj.(obj.oriList{i})), 2), 1);
            if ~isempty(idx)
                endIdx = min(idx(1)-1, endIdx);
            end
        end
    end
end