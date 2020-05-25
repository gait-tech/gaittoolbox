function out = calcSegAngAcc(obj, segs, frame)
	% Calculate the segment angular acceleration of grBody
	% 
	% Example: 
	%		out = obj.calcSegAngAcc({'LTIB', 'RTIB'})
	% 		out = struct('LTIB', n x 3, 'RTIB', n x 3)
	%
	% :param obj: this grBody
	% :param seg: cell array of segments to be calculated (if blank: all joints)
	% :param frame: out is expressed in this frame. B=body (default) or W=world
	% :return: out struct of velocities
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 18/02/20
    
    if nargin <= 1
        segs = obj.oriList;
    end
    if nargin <= 2
        frame = 'W';
    end
    validateattributes(frame,{'char'}, {'numel', 1})
    if ~(frame=='W')
        error(sprintf('Frame %s not supported\n', frame));
    end
    
    angvel = calcSegAngVel(obj, segs, frame);
    out = {}; fs = obj.fs;
    for i=1:length(segs)
        n = segs{i};
        angacc = diff(angvel.(n), 1, 1)*fs;
        out.(n) = [angacc(1,:); angacc];
    end
end