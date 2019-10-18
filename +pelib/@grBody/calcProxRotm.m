function prox = calcProxRotm(dist, angles, seq)
	% Calculate the proximal segment orientation rotation matrix
	% 
	% :param dist: quaternion orientation (n x 4) of the distal segment
	% :param angles: joint angles in seq order
	% :param seq: (default: YX'Z'')
	%
	% :return: dist quaternion orientation (n x 4) of the distal segment
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 3/20/19

    if nargin <= 2
        seq = 'YXZ';
    end
    
    relori = angle2quat(angles(:,1), angles(:,2), angles(:,3), seq);
    prox = quatmultiply(dist, quatconj(relori));
end