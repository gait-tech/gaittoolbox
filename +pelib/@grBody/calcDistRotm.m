function dist = calcDistRotm(prox, angles, seq)
	% Calculate the distal segment orientation rotation matrix
	% 
	% :param prox: quaternion orientation (n x 4) of the proximal segment
	% :param angles: joint angles in seq order
	% :param seq: (default: YX'Z'')
	%
	% :return: dist - quaternion orientation (n x 4) of the distal segment
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 3/20/19

    if nargin <= 2
        seq = 'YXZ';
    end
    
    relori = angle2quat(angles(:,1), angles(:,2), angles(:,3), seq);
    dist = quatmultiply(prox, relori);
end