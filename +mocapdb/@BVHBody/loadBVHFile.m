function obj = loadBVHFile(fname, unit)
	% Load BVH file and return an instance of BVHBody class
	%
	% :param fname: BVH file name
	% :param unit: data position unit (mm or inch)
	%
	% :return: instance of BVHBody class.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18

    if nargin > 1
        validStrings = ["mm", "inch"];
        unit = validatestring(unit, validStrings);
    else
        unit = "inch";
    end
    
    [bvh_data, time] = mocapdb.BVHBody.loadbvh(fname);
    
    obj = mocapdb.BVHBody('srcFileName', fname, 'posUnit', unit, ...
                         'frame', 'vicon');
        
    n = max(size(bvh_data));
    
    qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
    for i=1:n
        if length(bvh_data(i).name) >= 4 && strncmpi(bvh_data(i).name(end-3:end), 'base', 4)
            basename = bvh_data(i).name(1:end-4);
        else
            basename = bvh_data(i).name;
        end
        if isprop(obj, basename)
            v1 = basename;
            v2 = strcat('q', v1);
            
            if unit == "mm"
                obj.(v1) = (bvh_data(i).Dxyz*25.4)';
            else
                obj.(v1) = (bvh_data(i).Dxyz)';
            end
            obj.(v2) = rotm2quat(bvh_data(i).trans(1:3,1:3,:));
            obj.(v2) = quatmultiply(obj.(v2), qTCD2BM);
        end
    end
    
    obj.nSamples = length(bvh_data(1).Dxyz(1,:));
end