function out = getNonemptySeg(obj)
	% Get nonempty segment list
	%
	% :param obj: class XsensBody (self)
    %
	% :return: out - cell array of valid segments
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 Jun 9

    out = obj.segList;
    
    for i=length(obj.segList):-1:1
        n = out{i};
        if sum(size(obj.(n))) == 0
            out(i) = [];
        end
    end
end