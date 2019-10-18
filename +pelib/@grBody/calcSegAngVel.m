function out = calcSegAngVel(obj, segs, frame)
	% Calculate the segment angular velocity of grBody
	% 
	% Example: 
	%		out = obj.calcSegAngVel({'LTIB', 'RTIB'})
	% 		out = struct('LTIB', n x 3, 'RTIB', n x 3)
	%
	% :param obj: this grBody
	% :param seg: cell array of segments to be calculated (if blank: all joints)
	% :param frame: out is expressed in this frame. B=body (default) or W=world
	% :return: out struct of velocities
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 6/05/19

    if nargin <= 1
        segs = obj.oriList;
    end
    if nargin <= 2
        frame = 'B';
    end
    validateattributes(frame,{'char'}, {'numel', 1})
    if ~(frame=='B' || frame=='W')
        error(sprintf('Frame %s not supported\n', frame));
    end
    
    out = {}; fs = obj.fs;
    % Angular velocity of body in body frame as obtained from vicon input
    sIdx=1; eIdx=obj.nSamples;
    for i=1:length(segs)
        n = segs{i};
        w = quatmultiply(quatconj(obj.(n)(sIdx:eIdx, :)), ...
                                obj.(n)([sIdx+1:eIdx eIdx], :));
                            
        % this check is important as quat rep is a 2-1 mapping (i.e., q=-q)
        % I need angular velocity to be consistent (i.e., scalar part is
        % always positive)
        tmpIdx = w(:,1)<0;
        w(tmpIdx,:) = -w(tmpIdx,:);
        
        if frame == 'B'
            out.(n) = 2*w(:,2:4)*fs;
        elseif frame == 'W'
            out.(n) = quatrotate(quatconj(obj.(n)), 2*w(:,2:4)*fs);
        end
    end
end