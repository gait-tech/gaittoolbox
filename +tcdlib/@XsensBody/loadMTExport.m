% ======================================================================
%> @brief Load exported files from XSens MT manager (v4.8)
%>
%> Load exported files from XSens MT manager (v4.8)
%> Pelvis, L_UpLeg, R_UpLeg, L_LowLeg, R_LowLeg, L_Foot, R_Foot
%> 
%> options = struct('Pelvis', '00B40B91', ...
%> 'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%> 'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
%> 'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
%> Each field has a table with dimensions N x 13. The column of each row 
%> are quaternions[4], accelerometer [3], gyroscope [3], magnetometer [3]
%> 
%> For the TCD dataset, the accelerometer, gyroscope, magnetometer are in 
%> the sensor frame. quaternions tell the orientation relationship between 
%> sensor and world frame.
%>
%> @param name session name
%> @param options struct (body segment <-> sensor id)
%>
%> @retval data struct with the fields (dim N x 13) described. 
% ======================================================================
function obj = loadMTExport(name, options)
    % list of body segment names
    bs = fieldnames(options); 
    bsN = length(bs);
    
    obj = tcdlib.XsensBody('srcFileName', name, 'frame', 'sensor');
    nSamples = inf;
    for i = 1:bsN
        fpath = sprintf("%s-000_%s.txt", name, options.(bs{i}));
        T = readtable(fpath);
        obj.(bs{i}) = table([T.Quat_q0 T.Quat_q1 T.Quat_q2 T.Quat_q3], ...
                            [T.Acc_X T.Acc_Y T.Acc_Z], ...
                            [T.Gyr_X T.Gyr_Y T.Gyr_Z], ...
                            [T.Mag_X T.Mag_Y T.Mag_Z], ...
                           'VariableNames', {'ori', 'acc', 'gyr', 'mag'});
        Tsize = size(T);
        nSamples = min(nSamples, Tsize(1));
    end
    obj.nSamples = nSamples;
end