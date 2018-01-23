% ======================================================================
%> @brief Plot lower body segment length error
%> @author Luke Sy
%> 
%> @param body Body instance to be plotted
%> @param Hips , RFemur, LFemur, RTibia, LTibia
%>
%> @retval p plot object
% ======================================================================
function plotLowerBodySegmentLengthError(body, ...
    hips, lfemur, rfemur, ltibia, rtibia)
    
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