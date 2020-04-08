function out = calcCalibAnkleSensorW2PelvisWFromGyroSkewness(obj, DEGRANGE)
	% Calculate yaw offset from gyro skewness on sagital plane
	% 
	% :param obj: this XsensBody
	% :param calibS2B: XsensBody that transforms sensor frame to body frame
	%
	% :return: out - XsensBody class with adjustment sensor data
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 8/24/18

    if nargin <= 1, DEGRANGE = (0:1:359) - 180; end
    
    out = mocapdb.XsensBody();
    out.initializetoIdentity();
    
%     segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
%     ori = struct();
%     for i=1:length(segList)
%         n = segList{i};
%         ori.(n) = quatmultiply(obj.(n).ori, quatconj(calibS2B.(n).ori));
%     end
    
    segList = {'L_LowLeg', 'R_LowLeg'};   
    for i=1:length(segList)
        n = segList{i};
        
        err1 = zeros(length(DEGRANGE), 1);
%         err2 = zeros(length(DEGRANGE), 1);
%         isvalid = true(length(DEGRANGE), 1);
        
        for j=1:length(DEGRANGE)
            buf1 = quatmultiply(axang2quat([0 0 1 deg2rad(DEGRANGE(j))]), obj.(n).ori);
            buf2 = quatmultiply(quatconj(obj.Pelvis.ori), buf1);
            buf3 = quatrotate(quatconj(buf2), obj.(n).gyr);
            
            theta_y = buf3(:,2);
            err1(j) = skewness(theta_y);
%             err2(j) = peak2peak(theta_x);
            
%             theta_y0 = mean(theta_y(1:300));
%             [pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', mean(theta_y(1:300))+deg2rad(5));
%             [pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', mean(-theta_y(1:300))+deg2rad(5));
%             [pks1, locs1] = findpeaks(theta_y, 'MinPeakHeight', 0.75*theta_y0 + 0.25*max(theta_y));
%             [pks2, locs2] = findpeaks(-theta_y, 'MinPeakHeight', -0.75*theta_y0 + 0.25*max(-theta_y));
% 
%             if locs2(1) < locs1(1), isvalid(j) = false; end
        end
%         err1(~isvalid) = nan; 
%         err2(~isvalid) = nan;
        [~, I] = min(err1);
        out.(n).ori = axang2quat([0 0 1 deg2rad(DEGRANGE(I))]);
    end
end