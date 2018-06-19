% ======================================================================
%> @brief Transform BVHBody from vicon frame (default) to world frame
%>
%> @param obj this BVHBody
%> @param qR transformation quaternion (1 x 4) from vicon frame to world frame
%>
%> @retval out BVHBody in world frame.
% ======================================================================
function out = toWorldFrame(obj, qR)
    out = mocapdb.ViconBody();
    out.srcFileName = obj.srcFileName;
    out.frame = 'world';
    out.posUnit = obj.posUnit;
    out.fs = obj.fs;
    out.nSamples = obj.nSamples;
        
    posList = {'PELV', 'LFEP', 'LFEO', 'LTIO', 'LTOE', ...
        'RFEP', 'RFEO', 'RTIO', 'RTOE'};
           
    qR2 = quatconj(qR);
    for i=1:length(posList)
        out.(posList{i}) = quatrotate(qR2, obj.(posList{i}));
    end
    
    oriList = {'qRPV', 'qRTH', 'qLTH', 'qRSK', 'qLSK'};
	for i=1:length(oriList)
        out.(oriList{i}) = quatmultiply(qR, obj.(oriList{i}));
    end
end