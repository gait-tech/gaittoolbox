function out = generateBodyFromJointAngles(posMP, qOriMP, ...
    anglesLT, anglesRT, angleLK, angleRK, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, seq)
	% Generate grBody from pelvis pos, ori and joint angles
	%
	% :param posMP:  pelvis (root) position
	% :param qOriMP: pelvis (root) orientation
	% :param anglesLT: left  thigh/femur joint angles (X, Y, Z axis)
	% :param anglesRT: right thigh/femur joint angles (X, Y, Z axis)
	% :param angleLK: left  knee angle (X, Y, Z axis)
	% :param angleRK: right knee angle (X, Y, Z axis)
	% :param dPelvis: pelvis length
	% :param dLFemur: left  femur length
	% :param dRFemur: right femur length
	% :param dLTibia: left  tibia length
	% :param dRTibia: right tibia length
	% :param seq: (default: YX'Z'')
	%
	% :return: out new grBody in world frame
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    if nargin <= 11
        seq = 'YXZ';
    end
    seq2 = zeros(3, 1);
    for i=1:3, seq2(i) = seq(i)-'X'+1; end
    
    nSamples = size(posMP, 1);
    out = pelib.grBody();
    out.nSamples = size(posMP, 1);
    out.MIDPEL = posMP;
    out.qRPV = qOriMP;
    anglesLT = anglesLT .* [-1 -1 -1];
    anglesRT = anglesRT .* [1 -1 1];
    angleLK = angleLK .* [-1 1 -1];
    out.qLTH = pelib.grBody.calcDistRotm(out.qRPV, anglesLT(:, seq2), seq);
    out.qRTH = pelib.grBody.calcDistRotm(out.qRPV, anglesRT(:, seq2), seq);
    out.qLSK = pelib.grBody.calcDistRotm(out.qLTH, angleLK(:, seq2), seq);
    out.qRSK = pelib.grBody.calcDistRotm(out.qRTH, angleRK(:, seq2), seq);
    
    PELV_CS = quat2rotm(out.qRPV); PELV_Y = squeeze(PELV_CS(:,2,:))';
    out.LFEP = out.MIDPEL + dPelvis/2*PELV_Y;
    out.RFEP = out.MIDPEL - dPelvis/2*PELV_Y;
    LFEM_CS = quat2rotm(out.qLTH); LFEM_Z = squeeze(LFEM_CS(:,3,:))';
    RFEM_CS = quat2rotm(out.qRTH); RFEM_Z = squeeze(RFEM_CS(:,3,:))';
    out.LFEO = out.LFEP - dLFemur*LFEM_Z;
    out.RFEO = out.RFEP - dRFemur*RFEM_Z;
    LTIB_CS = quat2rotm(out.qLSK); LTIB_Z = squeeze(LTIB_CS(:,3,:))';
    RTIB_CS = quat2rotm(out.qRSK); RTIB_Z = squeeze(RTIB_CS(:,3,:))';
    out.LTIO = out.LFEO - dLTibia*LTIB_Z;
    out.RTIO = out.RFEO - dRTibia*RTIB_Z;
end