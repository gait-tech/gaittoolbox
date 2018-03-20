% ======================================================================
%> @brief change reference frame of grBody
%> Supported changes are vicon -> MIDPEL
%>
%> @param obj grBody (self)
%> @param ref reference frame
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = changeRefFrame(obj, ref)
    if nargin <= 1
        ref = 'MIDPEL0';
    end
    
    refMap = containers.Map({'MIDPEL'}, ...
        {'qRPV'});
    
    if strcmp(obj.frame, 'vicon') && isKey(refMap, ref)
        refPos = ref;
        refOri = refMap(ref);
    end
    
    out = copy(obj);
    out.frame = refPos;
    
    posList = obj.posList;
    for i=1:length(posList)
        out.(posList{i}) = quatrotate(obj.(refOri), ...
            obj.(posList{i})-obj.(refPos));
    end
    
    oriList = obj.oriList;
    for i=1:length(oriList)
        out.(oriList{i}) = quatmultiply(quatconj(obj.(refOri)), ...
            obj.(oriList{i}));
    end
end