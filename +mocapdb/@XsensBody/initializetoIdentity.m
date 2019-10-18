function initializetoIdentity(obj)
	% Initialize this XsensBody such that each rotation matrix is identity
	%
	% :param obj: this XsensBody
	%
	% .. Author: - Luke Sy (UNSW GSBME)
    
    for i=1:length(obj.segList)
        n = obj.segList{i};
        obj.(n).ori = rotm2quat(eye(3));
    end
end