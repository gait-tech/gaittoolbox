function out = calcCalibSB(obj, refBody, sIdx)
<<<<<<< HEAD
	% Calculate the calibration between sensor frame to body frame
=======
	% Calculate the calibration between sensor frame to body frame from one
	% frame (usually first frame).
	% 
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
	%
	% :param refBody: grBody class in world frame. (use data at index 1)
	% :param sIdx: index of obj to be used
	%
	% :return: out - Xsens with calibration data from sensor frame to body frame
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18

    if nargin <= 2
        sIdx = 1;
    end
    
    %% Variable initialization
    out = mocapdb.XsensBody('srcFileName', refBody.name, ...
                            'nSamples', 1, ...
                            'frame', 'calib');
    
    %% Calculation
<<<<<<< HEAD
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg', 'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val = {'qRPV', 'qLTH', 'qRTH', 'qLSK', 'qRSK', 'qLFT', 'qRFT'};
=======
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg', ...
           'L_LowLeg', 'R_LowLeg', 'L_LowLeg2', 'R_LowLeg2', ...
           'L_Foot', 'R_Foot'};
    val = {'qRPV', 'qLTH', 'qRTH', 'qLSK', 'qRSK', 'qLSK', 'qRSK', ...
           'qLFT', 'qRFT'};
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
    n = length(key);
    
    for i=1:n
        w_q_b = refBody.(val{i});
        if (~isempty(w_q_b) && ~isempty(obj.(key{i})))            
            w_q_s = obj.(key{i}).ori;
            out.(key{i}).ori = quatmultiply(quatconj(w_q_b(1,:)), w_q_s(sIdx,:));
        end
    end
end