function out = calcSegMeanRot(obj, seg, roty, idx)
	% Calculate the mean rotation of selected segments
	% 
	% Example 1:: 
    % 
	%	out = obj.calcSegMeanRot({'LTIO', 'RFEO'})
	%
	% :param obj: this object instance
    % :type obj: :class:`+pelib.@grBody`
	% :param seg: cell array string of body segments (see keys of :attr:`+pelib.@XsensBody.XsensBody.segList`)
    % :type seg: cell array
    % :param roty: radian rotation in y axis. Useful in adjusting between feet and other body segment orientation
    % :type roty: numeric in radians
	% :param idx: target indices. defaults to 1:obj.nSamples
    % :type idx: integer array
    %
	% :return: out - mean rotation
    % :rtype: length(idx) x 4
	%
	% .. Author: - Luke Sy (UNSW GSBME) Created 20/05/31
    if nargin <= 2
        roty = 0;
    end
    if nargin <= 3
        idx = 1:obj.nSamples;
    end
    if length(seg) <= 0
        error('Input segment must not be empty');
    end
    
    buf = quaternion(obj.(seg{1}).ori(idx,:));
    for i=2:length(seg)
        buf = [buf quaternion(obj.(seg{i}).ori(idx,:))];
    end
    qroty = quaternion([roty,0,0],'euler','YXZ','frame');
    out = compact(meanrot(buf, 2)*qroty);
end