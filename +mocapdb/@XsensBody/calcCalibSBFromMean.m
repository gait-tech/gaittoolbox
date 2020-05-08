function out = calcCalibSBFromMean(obj, refBody)
	% Calculate the calibration between sensor frame to body frame 
    % by taking the mean rotation offset between obj and refBody
	% 
    % Return B_q_S = argmin || ||
	% :param refBody: grBody class in world frame.
	%
	% :return: out - Xsens with calibration data from sensor frame to body frame
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 8
    
    %% Variable initialization
    out = mocapdb.XsensBody('srcFileName', refBody.name, ...
                            'nSamples', 1, ...
                            'frame', 'calib');
    
    %% Calculation
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg', 'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val = {'qRPV', 'qLTH', 'qRTH', 'qLSK', 'qRSK', 'qLFT', 'qRFT'};
    n = length(key);
    
    for i=1:n
        w_q_b = quaternion(refBody.(val{i}));
        if (~isempty(w_q_b) && ~isempty(obj.(key{i})))            
            w_q_s = quaternion(obj.(key{i}).ori(1:obj.nSamples,:));
            out.(key{i}).ori = compact(meanrot(w_q_b.conj().*w_q_s));
        end
    end
end