function out = changeRefFrame(obj, ref)
	% Change reference frame of grBody
	% Supported changes are vicon -> MIDPEL
	%
	% :param obj: grBody (self)
	% :param ref: reference frame
	%
	% :return: out - struct with the difference of pos and ori parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    if nargin <= 1
        ref = 'MIDPEL';
    end
    
    refMap = containers.Map({'MIDPEL'}, ...
        {'qRPV'});
    
    if strcmp(obj.frame, 'vicon') && isKey(refMap, ref)
        refPos = ref;
        refOri = refMap(ref);
    elseif strcmp(obj.frame, 'world') && isKey(refMap, ref)
        refPos = ref;
        refOri = refMap(ref);
    else
        error('The source frame or destination frame is not supported');
    end
    
    out = copy(obj);
    out.frame = refPos;
    
    posList = obj.posList;
    for i=1:length(posList)
        if (~isempty(obj.(posList{i})))
        out.(posList{i}) = quatrotate(obj.(refOri), ...
            obj.(posList{i})-obj.(refPos));
        end
    end
    
    oriList = obj.oriList;
    for i=1:length(oriList)
        if (~isempty(obj.(oriList{i})))
            out.(oriList{i}) = quatmultiply(quatconj(obj.(refOri)), ...
                obj.(oriList{i}));
        end
    end
end