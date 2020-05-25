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
        buf = diff(obj.(pts{i}), 2, 1)*fs*fs;
        out.(pts{i}) = [0 0 0; buf; buf(end,:)];
%         buf1 = obj.(pts{i});
%         buf2 = zeros(obj.nSamples,3);
%         for j=1:3
%             buf2(:,j) = gradient(gradient(buf1(:,j)));
%         end
%         out.(pts{i}) = buf2*fs*fs;
   end
end