function out = adjustFrame(obj, qR1, qR2)
	% adjust XsensBody frames by qR1 * obj.qB * qR2
	% 
	% :param obj: this XsensBody
	% :param qR1: quaternion 1
    % :param qR2: quaternion 2
	%
	% :return: out - XsensBody with adjusted convention.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 18

    out = obj.copy();
           
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n))
            out.(n).ori = quatmultiply(qR1, quatmultiply(obj.(n).ori, qR2));
        end
    end
end