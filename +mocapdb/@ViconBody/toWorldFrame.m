% ======================================================================
%> @brief Transform BVHBody from vicon frame (default) to world frame
%> @author Luke Sy (UNSW GSBME)
%> @date 24 Sept 2018
%>
%> @param obj this BVHBody
%> @param qR transformation quaternion (1 x 4) from vicon frame to world frame
%>
%> @retval out BVHBody in world frame.
% ======================================================================
function out = toWorldFrame(obj, qR)
    out = obj.copy();
    out.frame = 'world';
        
    posList = obj.posList;
           
    qR2 = quatconj(qR);
    for i=1:length(posList)
        out.(posList{i}) = quatrotate(qR2, obj.(posList{i}));
    end
    
    oriList = obj.oriList;
	for i=1:length(oriList)
        out.(oriList{i}) = quatmultiply(qR, obj.(oriList{i}));
    end
end