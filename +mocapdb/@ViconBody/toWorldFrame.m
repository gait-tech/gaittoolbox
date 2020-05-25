function out = toWorldFrame(obj, qR)
	% Transform BVHBody from vicon frame (default) to world frame
	%
	% :param obj: this ViconBody
	% :param qR: transformation quaternion (1 x 4) from vicon frame to world frame
	%
	% :return: out BVHBody in world frame.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/24/18

    out = obj.copy();
    out.frame = 'world';
        
    posList = obj.posList;
           
    qR2 = quatconj(qR);    
    posList = obj.posList;
    for i=1:length(posList)
        if(~isempty(obj.(posList{i})))
            out.(posList{i}) = quatrotate(qR2, obj.(posList{i}));
        end
    end
    
    oriList = obj.oriList;
    for i=1:length(oriList)
        if(~isempty(obj.(oriList{i})))
            out.(oriList{i}) = quatmultiply(qR, ...
                obj.(oriList{i}));
        end
    end
end