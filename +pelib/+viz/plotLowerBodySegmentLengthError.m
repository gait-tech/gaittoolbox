function plotLowerBodySegmentLengthError(body, ...
    hips, lfemur, rfemur, ltibia, rtibia)
    
	% Plot lower body segment length error
	% 
	% :param body: Body instance to be plotted
	% :param hips: left and right hip distance
	% :param lfemur: left femur length
	% :param rfemer: right femur length
	% :param ltibia: left tibia length
	% :param rtibia: right tibia length
	%
	% :return: p - plot object
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    errHips   = vecnorm(body.LFEP-body.RFEP, 2, 2) - hips;
    errLFemur = vecnorm(body.LFEP-body.LFEO, 2, 2) - lfemur;
    errRFemur = vecnorm(body.RFEP-body.RFEO, 2, 2) - rfemur;
    errLTibia = vecnorm(body.LFEO-body.LTIO, 2, 2) - ltibia;
    errRTibia = vecnorm(body.RFEO-body.RTIO, 2, 2) - rtibia;
    t = 1:length(errHips);
    
    p = plot(t, errHips, strcat(body.xyzColor{2},'-'), ...
             t, errLFemur, strcat(body.xyzColor{3},'--'), ...
             t, errRFemur, strcat(body.xyzColor{1},'--'), ...
             t, errLTibia, strcat(body.xyzColor{3},'-'), ...
             t, errRTibia, strcat(body.xyzColor{1},'-'));
    title('Segment Length Error');
    xlabel('Time'); ylabel(strcat('Error (', body.posUnit, ')'));
    legend('Hips', 'LFemur', 'RFemur', 'LTibia', 'RTibia');
end