% ======================================================================
%> @brief Load .sensors file from Total Capture dataset
%>
%> Load .sensors file and returns struct with fields:
%> Head, Sternum, Pelvis, L_UpArm, R_UpArm, L_LowArm, R_LowArm, L_UpLeg, 
%> R_UpLeg, L_LowLeg, R_LowLeg, L_Foot, R_Foot
%> 
%> Each field has a table with dimensions N x 13. The column of each row 
%> are quaternions[4], accelerometer [3], gyroscope [3], magnetometer [3]
%> 
%> For the TCD dataset, the accelerometer, gyroscope, magnetometer are in 
%> the sensor frame. quaternions tell the orientation relationship between 
%> sensor and world frame.
%>
%> @param fname .sensors filename
%>
%> @retval data struct with the fields (dim N x 13) described. 
% ======================================================================
function obj = loadSensorFile(fname)
    fileID = fopen(fname, 'r');
    buf = fscanf(fileID, "%d %d", 2);
    nSensors = buf(1);
    nFrames = buf(2);
    
    data0 = {};
    for i=1:nSensors
        data0{i} = [];
    end
    
    obj = mocapdb.XsensBody('srcFileName', fname, 'nSamples', nFrames, ...
                           'frame', 'sensor');
                           
    field = {};
    fscanf(fileID, "%d", 1);
    for i=1:nSensors
        buf = fscanf(fileID, "%s", 1);
        field{i} = buf;
        buf = fscanf(fileID, "%f", 13);
        data0{i}(1,:) = buf';
    end
       
    for i=2:nFrames
       fscanf(fileID, "%d", 1);
       for j=1:nSensors
           fscanf(fileID, "%s", 1);
           buf = fscanf(fileID, "%f", 13);
           data0{j}(i,:) = buf';
       end
    end
    fclose(fileID);
    
    for i=1:nSensors
        obj.(field{i}) = table(data0{i}(:,1:4), data0{i}(:,5:7), ...
                               data0{i}(:,8:10), data0{i}(:,11:13), ...
                               'VariableNames', ...
                               {'ori', 'acc', 'gyr', 'mag'});
    end 
end