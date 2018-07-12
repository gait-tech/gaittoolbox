% ======================================================================
%> @brief Load BVH file and return an instance of BVHBody class
%>
%> @param fname BVH file name
%> @param unit data position unit (mm or inch)
%>
%> @return instance of BVHBody class.
% ======================================================================
function obj = loadXsensBVHFile(fname, unit)
    [bvh_data, time] = mocapdb.BVHBody.loadbvh(fname);
    
    obj = mocapdb.BVHBody('srcFileName', fname, 'posUnit', unit, ...
                         'frame', 'world', 'fs', 1/(time(1)-time(2)));
    m = struct('Hips', 'Hips', 'Chest', 'Spine', 'Chest2', 'Spine1', ...
        'Chest3', 'Spine2', 'Chest4', 'Spine3', 'Neck', 'Neck', 'Head', 'Head', ...
        'RightCollar', 'RightShoulder', 'RightShoulder', 'RightArm', ...
        'RightElbow', 'RightForeArm', 'RightWrist', 'RightHand', ...
        'LeftCollar', 'LeftShoulder', 'LeftShoulder', 'LeftArm', ...
        'LeftElbow', 'LeftForeArm', 'LeftWrist', 'LeftHand', ...
        'RightHip', 'RightUpLeg', 'RightKnee', 'RightLeg',  ...
        'RightAnkle', 'RightFoot', 'RightToe', 'RightToe', ...
        'LeftHip', 'LeftUpLeg', 'LeftKnee', 'LeftLeg', ...
        'LeftAnkle', 'LeftFoot', 'LeftToe', 'LeftToe');
    n = max(size(bvh_data));
    
    qXsens2BM = rotm2quat([0 1 0; 0 0 1; 1 0 0]);
    for i=1:n
        if isfield(m, bvh_data(i).name)
            pname = m.(bvh_data(i).name);
        else
            pname = '';
        end
        if isprop(obj, pname)
            v1 = m.(bvh_data(i).name);
            v2 = strcat('q', v1);
            
            obj.(v1) = (bvh_data(i).Dxyz)';
            obj.(v2) = rotm2quat(bvh_data(i).trans(1:3,1:3,:));
%             obj.(v2) = rotm2quat(permute(bvh_data(i).trans(1:3,1:3,:), [2 1 3]));
            obj.(v2) = quatmultiply(obj.(v2), qXsens2BM);
        end
    end
    
    obj.nSamples = length(bvh_data(1).Dxyz(1,:));
end