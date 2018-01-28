% ======================================================================
%> @brief Load *_ori.txt and *_pos.txt file
%>
%> Load the orientations and positions obtained from the video.
%> Returns struct with fields similar to animatedata.
%>
%> @param fname .sensors filename
%>
%> @retval output struct with the fields described. 
% ======================================================================
function output = loadVideoOriPos(fname_ori, fname_pos)
    %% Check function input
    if nargin == 1
        fname = fname_ori;
        fname_ori = strcat(fname, '_ori.txt');
        fname_pos = strcat(fname, '_pos.txt');
    else
        validateattributes(fname_ori, 'string');
        validateattributes(fname_pos, 'string');
    end
    
    
    key = {'Hips', 'Spine', 'Neck', 'Head',...
        'RightShoulder', 'RightArm', 'RightForeArm', 'RightHand',...
        'LeftShoulder', 'LeftArm', 'LeftForeArm', 'LeftHand',...
        'RightUpLeg', 'RightLeg', 'RightFoot', 'RightToeBase',...
        'LeftUpLeg', 'LeftLeg', 'LeftFoot', 'LeftToeBase'};
    value = {'hips', 'spine', 'neck', 'head',...
        'rshoulder', 'rarm', 'rfarm', 'rhand',...
        'lshoulder', 'larm', 'lfarm', 'lhand',...
        'rupleg', 'rleg', 'rfoot', 'rtoe',...
        'lupleg', 'lleg', 'lfoot', 'ltoe'};
    m = containers.Map(key, value);
    output = struct();
    
    %% Load orientation data
    fileID_ori = fopen(fname_ori, 'r');
    buf = textscan(fileID_ori, '%s', 1, 'Delimiter', '\n');
    var_names = textscan(char(buf{1}), '%s');
    var_names = var_names{1};
    colN = length(var_names);
    data = fscanf(fileID_ori, '%f', [colN*4 Inf])';
    
    for i=1:colN
        if isKey(m, var_names{i})
            q = data(:,(i-1)*4+1:i*4);
            % rotm = quat2rotm(q);
%             rotm = q';
            output = setfield(output, strcat(m(var_names{i}), '_cs'), q(:,[4,1,2,3])');
        end
    end
    fclose(fileID_ori);
    
    %% Load position data
    fileID_pos = fopen(fname_pos, 'r');
    buf = textscan(fileID_pos, '%s', 1, 'Delimiter', '\n');
    var_names = textscan(char(buf{1}), '%s');
    var_names = var_names{1};
    colN = length(var_names);
    data = fscanf(fileID_pos, '%f', [colN*3 Inf]);
    
    for i=1:colN
        if isKey(m, var_names{i})
            pos = data((i-1)*3+1:i*3,:);
            output = setfield(output, m(var_names{i}), pos);
        end
    end
    
    fclose(fileID_pos);