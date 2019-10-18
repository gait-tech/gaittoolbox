function out = toViconFrame(obj, qR)
	% Transform XsensBody from world frame (default) to vicon frame
	% 
	% :param obj: this XsensBody
	% :param qR: XsensBody containing transformation quaternion (1 x 4) 
	%           from world frame to vicon frame
	%
	% :return: out - XsensBody in vicon frame.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/24/18

    out = obj.copy();
    out.frame = 'vicon';
           
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n))
            out.(n).ori = quatmultiply(qR.(n).ori, obj.(n).ori);
        end
    end
end