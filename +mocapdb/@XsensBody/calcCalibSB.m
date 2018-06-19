% ======================================================================
%> @brief Calculate the calibration between sensor frame to body frame
%>
%>
%> @param refBody grBody class in world frame
%>
%> @retval output Xsens with calibration data from sensor frame to body
% frame
% ======================================================================
function out = calcCalibSB(obj, refBody)
    %% Variable initialization
    out = mocapdb.XsensBody('srcFileName', refBody.name, ...
                            'nSamples', 1, ...
                            'frame', 'calib');
    
    %% Calculation
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg', 'L_LowLeg', 'R_LowLeg'};
    val = {'qRPV', 'qLTH', 'qRTH', 'qLSK', 'qRSK'};
    n = length(key);
    
    for i=1:n
        w_q_b = refBody.(val{i});
        w_q_s = obj.(key{i}).ori;
        out.(key{i}).ori = quatmultiply(quatconj(w_q_b(1,:)), w_q_s(1,:));
    end
end