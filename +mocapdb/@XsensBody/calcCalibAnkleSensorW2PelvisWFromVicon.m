% ======================================================================
%> @brief Calculate yaw offset from vicon data
%> @author Luke Sy (UNSW GSBME)
%> @date 06 Feb 2019
%>
%>
%> @param obj this XsensBody
%>
%> @return out XsensBody class with adjustment sensor data
% ======================================================================
function out = calcCalibAnkleSensorW2PelvisWFromVicon(obj, dataV)
    DEGRANGE = (0:0.1:359) - 180;
    out = mocapdb.XsensBody();
    
    segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
    
    nSamples = min(dataV.nSamples, obj.nSamples);
    sIdx = max(dataV.getStartIndex()+1, 100);
    eIdx = min(dataV.getEndIndex()-1, obj.nSamples-1);
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);

    viconBody = dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataV.fs, ...
                             'xyzColor', {'m', 'y', 'c'}});  
    refAcc = viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
    refAcc.MIDPEL = refAcc.MIDPEL(sIdx:eIdx,:);
    refAcc.LTIO = refAcc.LTIO(sIdx:eIdx,:);
    refAcc.RTIO = refAcc.RTIO(sIdx:eIdx,:);

    estAcc = {};
    estAcc.MIDPEL = quatrotate(quatconj(obj.Pelvis.ori), ...
                            obj.Pelvis.acc) - [0 0 9.81];
    estAcc.MIDPEL = estAcc.MIDPEL(sIdx:eIdx,:);
    estAcc.LTIO = quatrotate(quatconj(obj.L_LowLeg.ori), ...
                            obj.L_LowLeg.acc) - [0 0 9.81];
    estAcc.LTIO = estAcc.LTIO(sIdx:eIdx,:);
    estAcc.RTIO = quatrotate(quatconj(obj.R_LowLeg.ori), ...
                            obj.R_LowLeg.acc) - [0 0 9.81];
    estAcc.RTIO = estAcc.RTIO(sIdx:eIdx,:);

    [thetaP, errP] = findOptimalThetaBrute(refAcc.MIDPEL, estAcc.MIDPEL);
    [thetaL, errL] = findOptimalThetaBrute(refAcc.LTIO, estAcc.LTIO);
    [thetaR, errR] = findOptimalThetaBrute(refAcc.RTIO, estAcc.RTIO);

    [~, peakIdxP] = min(abs(DEGRANGE-thetaP)); 
    [~, peakIdxL] = min(abs(DEGRANGE-thetaL)); 
    [~, peakIdxR] = min(abs(DEGRANGE-thetaR)); 
    
    out.Pelvis.ori = axang2quat([0 0 1 deg2rad(-DEGRANGE(peakIdxP))]);
    out.L_LowLeg.ori = axang2quat([0 0 1 deg2rad(-DEGRANGE(peakIdxL))]);
    out.R_LowLeg.ori = axang2quat([0 0 1 deg2rad(-DEGRANGE(peakIdxR))]);
end

function [theta, err] = findOptimalThetaBrute(refAcc, estAcc)
    DEGRANGE = (0:0.1:359) - 180;
    err = zeros(length(DEGRANGE), 1);
    for i=1:length(DEGRANGE)
        estAcc2 = rotateThetaAlongZ(estAcc, DEGRANGE(i));
        err(i) = sqrt(mean(nanmean((refAcc(:,1:2) - estAcc2(:,1:2)).^2, 1)));
    end
    [M, I] = min(err);
    theta = DEGRANGE(I);
end

function newData = rotateThetaAlongZ(data, deg)
    R = [cosd(deg) sind(deg) 0;
         -sind(deg) cosd(deg) 0;
             0 0 1];
    newData = quatrotate(quatconj(rotm2quat(R)), data);
end