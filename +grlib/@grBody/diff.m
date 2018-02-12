% ======================================================================
%> @brief Returns the difference between grBody1 and grBody2
%>
%> @param obj1 grBody 1 (self)
%> @param obj1 grBody 2 (other)
%> @param seq orientation sequence
%>
%> @retval out struct with the difference of pos and ori parameters
% ======================================================================
function out = diff(obj1, obj2, seq)
    if nargin <= 2
        seq = 'YXZ';
    end
    
    out = struct;
    
    posList = obj1.posList;
    for i=1:length(posList)
        if length(obj1.(posList{i})) == 0 || length(obj2.(posList{i})) == 0
            out.(posList{i}) = [];
        else
            out.(posList{i}) = obj1.(posList{i}) - obj2.(posList{i});
        end
    end
    
    oriList = obj1.oriList;
    for i=1:length(oriList)
        if length(obj1.(oriList{i})) == 0 || length(obj2.(oriList{i})) == 0
            out.(oriList{i}) = [];
        else
            [r1 r2 r3] = quat2angle(obj1.(oriList{i}), seq);
            eul1 = [r1 r2 r3];
            [r1 r2 r3] = quat2angle(obj2.(oriList{i}), seq);
            eul2 = [r1 r2 r3];

            if obj1.oriUnit == 'deg'
                out.(oriList{i}) = grlib.anglediff(eul1, eul2)*180/pi;
            else
                out.(oriList{i}) = grlib.anglediff(eul1, eul2);
            end
        end
    end
end