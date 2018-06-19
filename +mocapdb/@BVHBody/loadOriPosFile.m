% ======================================================================
%> @brief Load BVH file and return an instance of BVHBody class
%>
%> @param fname_ori orientation file name
%> @param fname_pos position file name
%> @param unit data position unit (mm or inch)
%>
%> @return instance of BVHBody class.
% ======================================================================
function obj = loadOriPosFile(fname_ori, fname_pos, unit)
    %% Check function input
    if nargin == 1
        fname = fname_ori;
        fname_ori = strcat(fname, '_ori.txt');
        fname_pos = strcat(fname, '_pos.txt');
    else
        validateattributes(fname_ori, {'string', 'char'}, {});
        validateattributes(fname_pos, {'string', 'char'}, {});
    end
    
    if nargin > 2
        validStrings = ["mm", "inch"];
        unit = validatestring(unit, validStrings);
    else
        unit = "inch";
    end
    
    obj = mocapdb.BVHBody('srcFileName', {fname_ori, fname_pos}, ...
                         'posUnit', unit, 'frame', 'vicon');
    
    %% Load orientation data
    fileID_ori = fopen(fname_ori, 'r');
    buf = textscan(fileID_ori, '%s', 1, 'Delimiter', '\n');
    var_names = textscan(char(buf{1}), '%s');
    var_names = var_names{1};
    colN = length(var_names);
    data = fscanf(fileID_ori, '%f', [colN*4 Inf])';
    
    qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
    for i=1:colN
        if isprop(obj, var_names{i})
            v = strcat('q', var_names{i});
            q = data(:,(i-1)*4+1:i*4);
            % rotm = quat2rotm(q);
%             rotm = q';
            obj.(v) = quatmultiply(q(:,[4,1,2,3])', qTCD2BM);
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
        if isprop(obj, var_names{i})
            v = var_names{i};
            if unit == "mm"
                obj.(v) = (data((i-1)*3+1:i*3,:)*25.4)';
            else
                obj.(v) = (data((i-1)*3+1:i*3,:))';
            end
        end
    end
    
    obj.nSamples = length(obj.Hips(:,1));
    fclose(fileID_pos);
end