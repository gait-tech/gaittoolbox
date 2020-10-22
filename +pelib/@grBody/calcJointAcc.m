function out = calcJointAcc(obj, pts)
	% Calculate the joint acceleration of grBody
	% 
	% Example 1:: 
    % 
	%	out = obj.calcJointAcc({'LTIO', 'RFEO'})
	% 	% out = struct('LTIO', n x 3, 'RFEO', n x 3)
    %
    % Example 2:: 
    % 
    %   RFA_p_RF = randn(n, 4); RFA_p_RF(:,4) = 1;
	%	out = obj.calcJointVel(struct('LTIO': [0 0 1 1], 'RFT': RFA_p_RF))
	% 	% out = struct('LTIO', n x 3, 'RFT', n x 3)
    %
    % the values in the struct denote the position of target point
    % with respect the corresponding body frame (i.e., :math:`^{B}_{TP}p`).
    % For example, ``out.LTIO`` will contain the acceleration at point 
    % :math:`^{W}_{LS}\mathbf{T}` * ``pts.LTIO`` where :math:`^{W}_{LS}\mathbf{T}` 
    % is the transformation matrix of the left shanks wrt world frame and 
    % ``pts.LTIO`` is the position of the target point wrt the left 
    % shank frame, giving us the position of the target point wrt world 
    % frame. The same logic applies for the other keys (w/ the 
    % corresponding body frame) and values of input struct ``pts``.
	%
	% :param obj: this object instance
    % :type obj: :class:`+pelib.@grBody`
	% :param pts: cell array values or struct key indicate which body segment (see keys of :attr:`+pelib.@grBody.grBody.TMap`), struct value = 1 x 4 or n x 4 target point (3D homogenous point) wrt reference body, defaults to cell array with all joints
    % :type pts: cell array, struct, optional
	%
	% :return: out - struct of acceleration
    % :rtype: struct
	%
	% .. Author: - Luke Sy (UNSW GSBME) Created 18/09/24, Modified 20/05/31

    if nargin <= 1
        pts = obj.posList;
    end
    
    out = {}; fs = obj.fs;
    if iscell(pts)
        data = obj;
        keys = pts;
    else % struct
        data = obj.applyTMatrix(pts);
        keys = fieldnames(pts);
    end
    
    for i=1:length(keys)
        buf = diff(data.(keys{i}), 2, 1)*fs*fs;
        out.(keys{i}) = [0 0 0; buf; buf(end,:)];
%         buf1 = obj.(pts{i});
%         buf2 = zeros(obj.nSamples,3);
%         for j=1:3
%             buf2(:,j) = gradient(gradient(buf1(:,j)));
%         end
%         out.(pts{i}) = buf2*fs*fs;
    end
end