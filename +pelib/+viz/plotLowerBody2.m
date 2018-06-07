% ======================================================================
%> @brief Plot the lower body using the orientations (checking purposing)
%> @author Luke Sy
%>
%> Plot the lower body
%> 
%> @param body Body instance to be plotted
%> @param t time point to be plotted
%> @param showGround show ground
%> @retval p plot object
% ======================================================================
function p = plotLowerBody(body, t)    
    if nargin <= 2
        showGround = false;
    end
    dPelvis = norm(body.LFEP(t,:)-body.RFEP(t,:));
    dLFemur = norm(body.LFEP(t,:)-body.LFEO(t,:));
    dRFemur = norm(body.RFEP(t,:)-body.RFEO(t,:));
    dLTibia = norm(body.LFEO(t,:)-body.LTIO(t,:));
    dRTibia = norm(body.RFEO(t,:)-body.RTIO(t,:));
    
    MIDPEL = body.MIDPEL(t, :);
    LFEP = MIDPEL + 0.5*dPelvis*quatrotate(quatconj(body.qRPV(t, :)), [0 1 0]);
    RFEP = MIDPEL - 0.5*dPelvis*quatrotate(quatconj(body.qRPV(t, :)), [0 1 0]);
    LFEO = LFEP - dLFemur*quatrotate(quatconj(body.qLTH(t, :)), [0 0 1]);
    RFEO = RFEP - dRFemur*quatrotate(quatconj(body.qRTH(t, :)), [0 0 1]);
    LTIO = LFEO - dLTibia*quatrotate(quatconj(body.qLSK(t, :)), [0 0 1]);
    RTIO = RFEO - dRTibia*quatrotate(quatconj(body.qRSK(t, :)), [0 0 1]);
    
    pelvL = line([LFEP(1) MIDPEL(1) RFEP(1)], ...
                 [LFEP(2) MIDPEL(2) RFEP(2)], ...
                 [LFEP(3) MIDPEL(3) RFEP(3)], ...
                 'Color', char(body.xyzColor(2)), 'LineWidth', 2);
    pelvL.Marker = body.ptSymbol;
    lfemL = line([LFEP(1) LFEO(1)], ...
                 [LFEP(2) LFEO(2)], ...
                 [LFEP(3) LFEO(3)], ...
                 'Color', char(body.xyzColor(3)), 'LineWidth', 2, ...
                 'LineStyle', '--');
    rfemL = line([RFEP(1) RFEO(1)], ...
                 [RFEP(2) RFEO(2)], ...
                 [RFEP(3) RFEO(3)], ...
                 'Color', char(body.xyzColor(1)), 'LineWidth', 2, ...
                 'LineStyle', '--');
    ltibL = line([LFEO(1) LTIO(1)], ...
                 [LFEO(2) LTIO(2)], ...
                 [LFEO(3) LTIO(3)], ...
                 'Color', char(body.xyzColor(3)), 'LineWidth', 2);
    ltibL.Marker = body.ptSymbol;
    rtibL = line([RFEO(1) RTIO(1)], ...
                 [RFEO(2) RTIO(2)], ...
                 [RFEO(3) RTIO(3)],...
                 'Color', char(body.xyzColor(1)), 'LineWidth', 2);
    rtibL.Marker = body.ptSymbol;
    if showGround
        [ptsX, ptsY, ptsZ] = body.groundCoordinates();
        patch(ptsX, ptsY, ptsZ, 'k', 'FaceAlpha', .5, 'LineStyle', ':')
    end
end