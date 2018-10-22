% ======================================================================
%> @brief Calculate rotation matrix of sensor world frame with respect vicon frame
%> @author Luke Sy (UNSW GSBME)
%> @date 24 Sept 2018
%>
%> @param viconFName vicon calib mat file name
%> @param xsensFName xsens MT export file name
%> @param options XsensBody load configuration (see loadMTExport)
%> @param idx [OPTIONAL] index of xsens ori to be used in calculation
%>            (default = 1)
%>
%> @return obj XsensBody class qOri world to vicon frame of each sensor ({}^V_W q)
% ======================================================================
function obj = loadCalibSensorW2V(viconFName, xsensFName, options, idx)
    if nargin <= 3, idx = 5; end
    
    obj = mocapdb.XsensBody();
    
    w_s__q = mocapdb.XsensBody.loadMTExport(xsensFName, options);
    load(viconFName);
    
    m = struct('Pelvis', 'RPV', 'L_UpLeg', 'LTH', 'R_UpLeg', 'RTH', ...
            'L_LowLeg', 'LAK', 'R_LowLeg', 'RAK', ...
            'L_Foot', 'LFT', 'R_Foot', 'RFT');
    
    for i=1:length(w_s__q.segList)
        n = w_s__q.segList(i); n = n{1};
        
        if isfield(m, n)
            pN = eval(sprintf('%s_SensN', m.(n)));
            pS = eval(sprintf('%s_SensS', m.(n)));
            pW = eval(sprintf('%s_SensW', m.(n)));

%             x = nanmean(pN, 1) -  nanmean(pS, 1);
            x = pN(idx, :) -  pS(idx, :);
            x = x / norm(x);

%             y = nanmean(pW, 1) - nanmean(pN, 1);
            y = pW(idx, :) - pN(idx, :);
            y = y - dot(x, y)*x;
            y = y / norm(y);

            z = cross(x, y);
            z = z / norm(z);

            w_v__q = rotm2quat([x' y' z']);
            obj.(n).ori = quatmultiply(w_v__q, quatconj(w_s__q.(n).ori(idx, :)));
        end
    end
end