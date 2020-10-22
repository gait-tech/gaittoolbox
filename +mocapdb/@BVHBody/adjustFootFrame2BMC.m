function out = adjustFootFrame2BMC(obj)
	% Adjust the foot frame from animation convention to biomechanical
    % convention. See
    % `link <https://docs.vicon.com/display/Nexus29/Lower+body+kinematics#Lowerbodykinematics-Foot>`_
	% for description of convention.
    %
	% :param obj: this object instance
    % :type obj: :class:`+mocapdb.@BVHBody`
	%
	% :return: out - new object instance with correct foot
    % :rtype: :class:`+mocapdb.@BVHBody`
	%
	% .. Author: - Luke Sy (UNSW GSBME) Created Modified 20/06/01
    out = obj.copy();
    
    out.qLeftFoot = quatmultiply(out.qLeftFoot, axang2quat([0 1 0 pi/2]) );
    out.qRightFoot = quatmultiply(out.qRightFoot, axang2quat([0 1 0 pi/2]) );
    
    LFT_len = dot(quatrotate(quatconj(out.qLeftFoot(1,:)), [0 0 1]), ...
                  out.LeftToe(1,:)-out.LeftFoot(1,:));
    RFT_len = dot(quatrotate(quatconj(out.qRightFoot(1,:)), [0 0 1]), ...
                  out.RightToe(1,:)-out.RightFoot(1,:));
    LFT_z = quatrotate(quatconj(out.qLeftFoot), [0 0 LFT_len]);
    RFT_z = quatrotate(quatconj(out.qRightFoot), [0 0 RFT_len]);
    out.LeftToe = out.LeftFoot + LFT_z;
    out.RightToe = out.RightFoot + RFT_z;
end