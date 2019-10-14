function out = calcMPLARAdist(obj)    
	% Calculate the midpelvis, left ankle, and right ankle distances of grBody
	% 
	% Example: 
	% 		out = obj.calcMPLARAdist(100)
	% 		out = struct('MPLA', n x 1, 'MPRA', n x 1, 'LARA', n x 1)
	%
	% :param obj: this grBody
	%
	% :return: out struct of distances
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 3/15/19

    out = {};
    out.MPLA = vecnorm(obj.MIDPEL-obj.LTIO, 2, 2);
    out.MPRA = vecnorm(obj.MIDPEL-obj.RTIO, 2, 2);
    out.LARA = vecnorm(obj.LTIO-obj.RTIO, 2, 2);
end