% ======================================================================
%> @brief Transform BVHBody from vicon frame (default) to world frame
%>
%> @param obj this BVHBody
%> @param qR transformation quaternion (1 x 4) from vicon frame to world frame
%>
%> @retval out BVHBody in world frame.
% ======================================================================
function out = toWorldFrame(obj, qR)
    out = tcdlib.BVHBody();
    out.srcFileName = obj.srcFileName;
    out.frame = 'world';
    out.posUnit = obj.posUnit;
    out.nSamples = obj.nSamples;
        
    posList = {'Hips', 'Spine', 'Spine1', 'Spine2', 'Spine3', 'Neck', 'Head', ...
        'RightShoulder', 'RightArm', 'RightForeArm', 'RightHand', ...
        'LeftShoulder', 'LeftArm', 'LeftForeArm', 'LeftHand', ...
        'RightUpLeg', 'RightLeg', 'RightFoot', ...
        'LeftUpLeg', 'LeftLeg', 'LeftFoot'};
           
    qR2 = quatconj(qR);
    for i=1:length(posList)
        out.(posList{i}) = quatrotate(qR2, obj.(posList{i}));
    end
    
    oriList = {'qHips', 'qSpine', 'qSpine1', 'qSpine2', 'qSpine3', ...
               'qNeck', 'qHead', ...
               'qRightShoulder', 'qRightArm', 'qRightForeArm', 'qRightHand', ...
               'qLeftShoulder', 'qLeftArm', 'qLeftForeArm', 'qLeftHand', ...
               'qRightUpLeg', 'qRightLeg', 'qRightFoot', ...
               'qLeftUpLeg', 'qLeftLeg', 'qLeftFoot'};
	for i=1:length(oriList)
        out.(oriList{i}) = quatmultiply(qR, obj.(oriList{i}));
    end
end