function out = toWorldFrame(obj, qR)
	% Transform BVHBody from vicon frame (default) to world frame
	%
	% :param obj: this BVHBody
	% :param qR: transformation quaternion (1 x 4) from vicon frame to world frame
	%
	% :return: out - BVHBody in world frame.
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    out = mocapdb.BVHBody();
    out.srcFileName = obj.srcFileName;
    out.frame = 'world';
    out.posUnit = obj.posUnit;
    out.nSamples = obj.nSamples;
        
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