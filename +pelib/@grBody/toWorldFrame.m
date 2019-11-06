function out = toWorldFrame(obj, pos, ori)
	% Change reference frame of grBody from MIDPEL frame to world frame
	% Supported changes are vicon -> MIDPEL
	%
	% :param obj: grBody (self)
	% :param pos: pelvis (root) position
	% :param ori: pelvis (root) orientation
	%
	% :return: out - new grBody in world frame
	%
	% .. Author: - Luke Sy (UNSW GSBME)
    
    refMap = containers.Map({'MIDPEL'}, ...
        {'qRPV'});
    
    if ~strcmp(obj.frame, 'MIDPEL') && isKey(refMap, ref)
        error('The source frame is not supported');    
    end
    
    out = copy(obj);
    out.frame = 'world';
    
    posList = obj.posList;
    for i=1:length(posList)
        if(~isempty(obj.(posList{i})))
            out.(posList{i}) = quatrotate(quatconj(ori), ...
                obj.(posList{i})) + pos;
        end
    end
    
    oriList = obj.oriList;
    for i=1:length(oriList)
        if(~isempty(obj.(oriList{i})))
            out.(oriList{i}) = quatmultiply(ori, ...
                obj.(oriList{i}));
        end
    end
end