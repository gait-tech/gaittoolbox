% ======================================================================
%> @brief Calculate knee angle from 3P dist
%> @author Luke Sy (UNSW GSBME)
%> @date 22 Mar 2019
%>
%> @param prox quaternion orientation (n x 4) of the proximal segment
%> @param angles joint angles in seq order
%> @param seq (default: YX'Z'')
%>
%> @retval dist quaternion orientation (n x 4) of the distal segment
% ======================================================================
function [alphaLK, alphaRK, d] = calcKneeAnglesFromMPLARADist(...
                            PELV_CS, LTIB_CS, RTIB_CS, ...
                            dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, ...
                            dMPLADist, dMPRADist)
    dMPLADist = min(max(dMPLADist, 0), ...
        norm(dPelvis/2*PELV_CS(:,2)+(dLFemur+dLTibia)*LTIB_CS(:,3)));
    dMPRADist = min(max(dMPRADist, 0), ...
        norm(-dPelvis/2*PELV_CS(:,2)+(dRFemur+dRTibia)*RTIB_CS(:,3)));
    
    % calculate temp vector wL and wR
    wL =  dPelvis/2*PELV_CS(:,2) - dLTibia*LTIB_CS(:,3);
    wR = -dPelvis/2*PELV_CS(:,2) - dRTibia*RTIB_CS(:,3);
    
    aL = -2*dLFemur*dot(wL,LTIB_CS(:,3));
    bL =  2*dLFemur*dot(wL,LTIB_CS(:,1));
    cL = dMPLADist.^2 - dot(wL, wL) - dLFemur.^2;
    aR = -2*dRFemur*dot(wR,RTIB_CS(:,3));
    bR =  2*dRFemur*dot(wR,RTIB_CS(:,1));
    cR = dMPRADist.^2 - dot(wR, wR) - dRFemur.^2;
    
    % calculate alphaLK and alphaRK
    v1 = max(min(1,abs((aL.*cL+bL.*sqrt(aL.^2+bL.^2-cL.^2))/(aL.^2+bL.^2))),0);
    v2 = max(min(1,abs((aL.*cL-bL.*sqrt(aL.^2+bL.^2-cL.^2))/(aL.^2+bL.^2))),0);
    v3 = max(min(1,abs((aR.*cR+bR.*sqrt(aR.^2+bR.^2-cR.^2))/(aR.^2+bR.^2))),0);
    v4 = max(min(1,abs((aR.*cR-bR.*sqrt(aR.^2+bR.^2-cR.^2))/(aR.^2+bR.^2))),0);
    alphaLK = [ acos(v1), acos(v2) ];
    alphaRK = [ acos(v3), acos(v4) ];

    % test and return that value that satisfies the equation the most
%     diffDistLK = abs(vecnorm(dPelvis/2*PELV_CS(:,2) ...
%                          -dLFemur*cos(alphaLK).*LTIB_CS(:,3) ...
%                          +dLFemur*sin(alphaLK).*LTIB_CS(:,1) ...
%                          -dLTibia*LTIB_CS(:,3), 2, 1) - dMPLADist);
%     diffDistRK = abs(vecnorm(-dPelvis/2*PELV_CS(:,2)+ ...
%                          -dRFemur*cos(alphaRK).*RTIB_CS(:,3) ...
%                          +dRFemur*sin(alphaRK).*RTIB_CS(:,1) ...
%                          -dRTibia*RTIB_CS(:,3), 2, 1) - dMPRADist);
%     [~, iLK] = min(diffDistLK);
%     [~, iRK] = min(diffDistRK);
%     alphaLK = alphaLK(iLK); alphaRK = alphaRK(iRK);
end