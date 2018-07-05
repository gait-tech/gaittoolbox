% ======================================================================
%> @brief Load exported files from XSens MVN Studio 4.4 (.mvnx)
%>
%> Load exported files from XSens MT manager (v4.8)
%> Pelvis, L_UpLeg, R_UpLeg, L_LowLeg, R_LowLeg, L_Foot, R_Foot
%> 
%> options = struct('Pelvis', 'Pelvis', ...
%> 'L_UpLeg', 'LeftUpperLeg', 'R_UpLeg', 'RightUpperLeg', ...
%> 'L_LowLeg', 'prop', 'R_LowLeg', 'prop_1', ...
%> 'L_Foot', 'LeftFoot', 'R_Foot', 'RightFoot');
%> Each field has a table with dimensions N x 13. The column of each row 
%> are quaternions[4], accelerometer [3], gyroscope [3], magnetometer [3]
%> 
%> Accelerometer, gyroscope, magnetometer are in the sensor frame. 
%> Quaternions tell the orientation relationship between sensor and world frame.
%>
%> @param fname file name
%> @param options struct (body segment <-> sensor id)
%>
%> @retval XsensBody
% ======================================================================
function obj = loadMVNX(fname, options)
    if nargin <= 1
        options = struct('Pelvis', 'Pelvis', ...
            'L_UpLeg', 'LeftUpperLeg', 'R_UpLeg', 'RightUpperLeg', ...
            'L_LowLeg', 'LeftAnkle', 'R_LowLeg', 'RightAnkle', ...
            'L_Foot', 'LeftFoot', 'R_Foot', 'RightFoot');
    end

    % load data
    tree = mocapdb.XsensBody.load_mvnx(fname);
    nSamples = tree.subject.frames.frame(end).index+1;
    obj = mocapdb.XsensBody('srcFileName', fname, ...
                            'fs', tree.subject.frameRate, ...
                            'nSamples', nSamples, 'frame', 'sensor');

    %creates a struct with sensor data
    if isfield(tree.subject,'sensors') && isstruct(tree.subject.sensors)
        sensorData = tree.subject.sensors.sensor;
    else
        error('tree.subject.sensors does not exist');
    end
    
    % build idx which indicates which columns from mvnx should be exported
    sensorData2 = matlab.lang.makeUniqueStrings({sensorData.label});
    M = containers.Map(sensorData2, 1:length(sensorData2));
    optionFieldNames = fieldnames(options);
    nSensors = length(optionFieldNames);
    idx3 = []; idx4 = [];
    for i=1:nSensors
        j = M(options.(optionFieldNames{i}));
        idx3(end+1:end+3) = (j-1)*3+1:j*3;
        idx4(end+1:end+4) = (j-1)*4+1:j*4;
        remove(M, options.(optionFieldNames{i}));
    end
    
    acc = zeros(nSamples, nSensors*3);
    gyr = zeros(nSamples, nSensors*3);
    mag = zeros(nSamples, nSensors*3);
    ori = zeros(nSamples, nSensors*4);
    i = 1; j = 1;
    while i < length(tree.subject.frames.frame)
        if tree.subject.frames.frame(i).index
            acc(j, :) = tree.subject.frames.frame(i).sensorAcceleration(idx3);
            gyr(j, :) = tree.subject.frames.frame(i).sensorAngularVelocity(idx3);
            mag(j, :) = tree.subject.frames.frame(i).sensorMagneticField(idx3);
            ori(j, :) = tree.subject.frames.frame(i).sensorOrientation(idx4);
            j = j + 1;
        end
        i = i + 1;
    end

    for i=1:nSensors
        tIdx3 = (i-1)*3+1:i*3;
        tIdx4 = (i-1)*4+1:i*4;
        obj.(optionFieldNames{i}) = table(ori(:,tIdx4), acc(:,tIdx3), ...
            gyr(:,tIdx3), mag(:,tIdx3), 'VariableNames', ...
            {'ori', 'acc', 'gyr', 'mag'});
    end 
end
