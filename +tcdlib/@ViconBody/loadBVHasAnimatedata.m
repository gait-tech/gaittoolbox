function data = loadBVHasAnimatedata(fname, unit)
    if nargin > 1
        validStrings = ["mm", "inch"];
        unit = validatestring(unit, validStrings);
    else
        unit = "inch";
    end
    
    [bvh_data, time] = loadbvh(fname);
    
    data = struct('ltoe', [], 'lfoot', [], 'lleg', [], 'lupleg', [],...
        'rtoe', [], 'rfoot', [], 'rleg', [], 'rupleg', [],...
        'hips', [], 'spine', [], 'neck', [], 'head', [],...
        'rshoulder', [], 'rarm', [], 'rfarm', [], 'rhand', [],...
        'lshoulder', [], 'larm', [], 'lfarm', [], 'lhand', []);
    
    n = max(size(bvh_data));
    
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
    
    for i=1:n
        if isKey(m, bvh_data(i).name)
            v = m(bvh_data(i).name);
            
            if unit == "mm"
                data = setfield(data, v, bvh_data(i).Dxyz*25.4);
            else
                data = setfield(data, v, bvh_data(i).Dxyz);
            end
            data = setfield(data, strcat(v, '_cs'),...
                            bvh_data(i).trans(1:3,1:3,:));
        end
    end