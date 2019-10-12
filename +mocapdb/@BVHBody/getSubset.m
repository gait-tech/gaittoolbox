function out = getSubset(obj, idx)
	% Get subset of bvh body
	% 
	% Example:
	%      out = obj.getSubset(5:100);
	%
	% :param obj: class BVHBody (self)
	% :param idx: indices of data to be included in out
	% :return: out - BVHBody class whose data only includes the rows in idx
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 11/14/18

    validateattributes(idx, {'numeric'}, {});
    
    out = obj.copy();
    segList = [obj.posList obj.oriList];
    
    for i=1:length(segList)
        n = segList{i};
        if sum(size(out.(n))) ~= 0
            out.(n) = out.(n)(idx, :);
        end
    end
    out.nSamples = length(idx);
end