function out = calcCalibAnkleSensorW2PelvisWFromVicon(obj, dataV)
	% Calculate yaw offset from vicon data
	% 
	% :param obj: this XsensBody
	%
	% :return: out XsensBody class with adjustment sensor data
	% 
	% .. Author: - Luke Sy (UNSW GSBME) - 02/06/19

    DEGRANGE = (0:0.1:359) - 180;
    out = mocapdb.XsensBody();
    
    segList = ["Pelvis", "L_LowLeg", "R_LowLeg", "L_Foot", "R_Foot"];
    
    nSamples = min(dataV.nSamples, obj.nSamples);
    sIdx = max(dataV.getStartIndex()+1, 100);
    eIdx = min(dataV.getEndIndex()-1, obj.nSamples-1);
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);

    viconBody = dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataV.fs, ...
                             'xyzColor', {'m', 'y', 'c'}});  
    refAcc = viconBody.calcJointAcc(struct('MIDPEL', [0 0 0 1], ...
        'LTIO', [0 0 0 1], 'RTIO', [0 0 0 1], ...
        'LFT', [0 0 0.5*viconBody.calcLFootLength(sIdx) 1], ...
        'RFT', [0 0 0.5*viconBody.calcRFootLength(sIdx) 1] ));
    refAcc.Pelvis = refAcc.MIDPEL(sIdx:eIdx,:);
    refAcc.L_LowLeg = refAcc.LTIO(sIdx:eIdx,:);
    refAcc.R_LowLeg = refAcc.RTIO(sIdx:eIdx,:);
    refAcc.L_Foot = refAcc.LFT(sIdx:eIdx,:);
    refAcc.R_Foot = refAcc.RFT(sIdx:eIdx,:);

    estAcc = {};
    for i = segList
        estAcc.(i) = quatrotate(quatconj(obj.(i).ori), ...
                            obj.(i).acc) - [0 0 9.81];
        estAcc.(i) = estAcc.(i)(sIdx:eIdx,:);

        [theta, err] = findOptimalThetaBrute(refAcc.(i), estAcc.(i));
        [~, peakIdx] = min(abs(DEGRANGE-theta)); 
        out.(i).ori = axang2quat([0 0 1 deg2rad(-DEGRANGE(peakIdx))]);
    end
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