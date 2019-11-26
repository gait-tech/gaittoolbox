function out = calcJointAcc(obj, pts)
	% Calculate the joint acceleration of grBody
	% 
	% Example: 
	% 		out = obj.calcJointVel(100, {'LTIO', 'RFEO'})
	% 		out = struct('LTIO', n x 3, 'RFEO', n x 3)
	%
	% :param obj: this grBody
	% :param pts: cell array of joints to be calculated (if blank: all joints)
	%
	% :return: out struct of velocities
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/24/18

    if nargin <= 1
        pts = obj.posList;
    end
    
    out = {}; fs = obj.fs;
    for i=1:length(pts)
        n = pts{i};
        out.(n) = [0 0 0; diff(obj.(n), 2, 1)*fs*fs; 0 0 0];
    end
end