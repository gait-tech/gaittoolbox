% ======================================================================
%> @brief Calculate yaw offset from high RoM on sagital plane
%> @author Luke Sy (UNSW GSBME)
%> @date 22 Oct 2018
%>
%> @param obj this XsensBody
%> @param calibS2B XsensBody that transforms sensor frame to body frame
%>
%> @return out XsensBody class with adjustment sensor data
% ======================================================================
function out = calcCalibAnkleSensorW2PelvisWFromROM(obj, calibS2B, DEGRANGE)
    if nargin <= 2, DEGRANGE = (0:1:359) - 180; end
    
    out = mocapdb.XsensBody();
    out.initializetoIdentity();
    
    segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
    ori = struct();
    for i=1:length(segList)
        n = segList{i};
        ori.(n) = quatmultiply(obj.(n).ori, quatconj(calibS2B.(n).ori));
    end
    
    segListPV = {'L_LowLeg', 'R_LowLeg'};   
    for i=1:length(segListPV)
        n = segListPV{i};
        
        err1 = zeros(length(DEGRANGE), 1);
%         err2 = zeros(length(DEGRANGE), 1);
        isvalid = true(length(DEGRANGE), 1);
        
        warning('off', 'signal:findpeaks:largeMinPeakHeight');
        for j=1:length(DEGRANGE)
            buf1 = quatmultiply(axang2quat([0 0 1 deg2rad(DEGRANGE(j))]), ori.(n));
            buf2 = quatmultiply(quatconj(ori.Pelvis), buf1);
            [theta_y theta_x theta_z] = quat2angle(buf2, 'YXZ');
%             err(j) = sum(theta_y.^2);
            err1(j) = peak2peak(theta_y);
%             err2(j) = peak2peak(theta_x);
            
            theta_y0 = mean(theta_y(1:300));
            [pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', mean(theta_y(1:300))+deg2rad(10));
            [pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', mean(-theta_y(1:300))+deg2rad(10));
%             [pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', 0.1*theta_y0 + 0.1*max(theta_y));
%             [pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', -0.1*theta_y0 + 0.1*max(-theta_y));
            
            if length(pks1) == 0 || length(pks2) == 0
                isvalid(j) = false;
            else
                isvalid(j) = locs2(1) >= locs1(1);
            end
        end
        warning('on', 'signal:findpeaks:largeMinPeakHeight');
        
        [pks, locs] = findpeaks(err1);
        locs = locs(isvalid(locs));
        out.(n).ori = axang2quat([0 0 1 deg2rad(DEGRANGE(locs))]);
    end
end