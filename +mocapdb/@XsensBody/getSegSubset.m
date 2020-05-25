function out = getSegSubset(obj, segList)
	% Get segment subset of xsens measurements
	%
	% Example:
    %      segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'}
	%      out = obj.getSegSubset(segList);
	%
	% :param obj: class XsensBody (self)
	% :param segList: list of segments
	% :return: out - XsensBody class whose data only includes the rows in idx
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/20/18
    
    out = obj.copy();
    for i=1:length(segList)
        out.(segList{i}) = obj.(segList{i});
    end
end