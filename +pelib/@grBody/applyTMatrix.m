function out = applyTMatrix(obj, pts)
	% Apply the transformation matrix of grBody to points
	% 
    % Example 1:: 
    % 
	%	out = obj.applyTMatrix(struct('LTIO': [0 0 1 1], 'RFEO': [1 0 0 1]))
	% 	% out = struct('LTIO', n x 3, 'RFEO', n x 3)
    %
    % ``out.LTIO`` will contain :math:`^{W}_{LS}\mathbf{T}` * ``pts.LTIO``. 
    % The same logic applies for the other values (w/ the corresponding 
    % body frame) of input cell array ``pts``.
	%
	% :param obj: this object instance
    % :type obj: :class:`+pelib.@grBody`
	% :param pts: cell array values indicate which body segment (see keys of :attr:`+pelib.@grBody.grBody.TMap`), struct value = 1 x 4 or n x 4 target point (3D homogenous point) wrt reference body
    % :type pts: cell array
	%
	% :return: out - struct of points (n x 3)
    % :rtype: struct
	%
	% .. Author: - Luke Sy (UNSW GSBME) Created 20/05/31

    if nargin <= 1
        pts = obj.posList;
    end
    
    out = struct();
    keys = fieldnames(pts);
    data = obj.getTMatrix(keys);

    for i=1:length(keys)
        n = keys{i};
        bufOut = zeros(obj.nSamples, 4);
        if size(pts.(n), 1) == 1
            bufPts = pts.(n)';
            for j=1:obj.nSamples
                bufOut(j,:) = (data.(n)(:,:,j)*bufPts)';
            end
        elseif size(pts.(n), 1) == obj.nSamples
            bufPts = pts.(n)';
            for j=1:obj.nSamples
                bufOut(j,:) = (data.(n)(:,:,j)*bufPts(:,j))';
            end
        else
            error("Invalid size of pts.%s = %d", n, size(pts.(n), 2));
        end
        out.(n) = bufOut(:,1:3);
    end        
end