function out = getTMatrix(obj, pts, idx)
	% Calculate the transformation matrix of grBody
	% 
	% Example 1:: 
    % 
	%	out = obj.getTMatrix({'LTIO', 'RFEO'})
	% 	% out = struct('LTIO', 4 x 4 x n, 'RFEO', 4 x 4 x n)
    %
    % ``out.LTIO`` will contain the transformation matrix % :math:`^{W}_{LS}\mathbf{T}`
    % which containst that orientation of the left shank centered at ``pts.LTIO``
    % The same logic applies for the other values (w/ the corresponding 
    % body frame) of input cell array ``pts``.
	%
    % Example 2:: 
    % 
	%	out = obj.getTMatrix('LTIO', 1:3)
	% 	
	% :param obj: this object instance
    % :type obj: :class:`+pelib.@grBody`
	% :param pts: cell array values or characters indicating which body segment (see keys of :attr:`+pelib.@grBody.grBody.TMap`), defaults to cell array with all joints
    % :type pts: cell array, characters, optional
    % :param idx: target indices. defaults to 1:obj.nSamples
    % :type pts: integer array
	%
	% :return: out - struct of transformation matrices (4 x 4 x n) or just the transformation matrix
    % :rtype: struct, matrix
	%
	% .. Author: - Luke Sy (UNSW GSBME) Created 20/05/31

    if nargin <= 1
        pts = obj.posList;
    end
    if nargin <= 2
        idx = 1:obj.nSamples;
    end
    map = obj.TMap;
    
    if iscell(pts)
        out = struct();

        for i=1:length(pts)
            n = pts{i};
            buf = zeros(4,4,length(idx));
            buf(1:3,1:3,:) = quat2rotm(obj.(map.(n).ori)(idx,:));
            buf(1:3,4,:) = obj.(map.(n).trans)(idx,:)';
            buf(4,4,:) = 1;
            out.(n) = buf;
        end
    elseif ischar(pts) || isstring(pts)
        out = zeros(4,4,length(idx));
        out(1:3,1:3,:) = quat2rotm(obj.(map.(pts).ori)(idx,:));
        out(1:3,4,:) = obj.(map.(pts).trans)(idx,:)';
        out(4,4,:) = 1;
    else
        error('Input type %s for pts is not supported', class(pts));
    end
end