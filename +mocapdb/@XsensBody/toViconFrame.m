% ======================================================================
%> @brief Transform XsensBody from world frame (default) to vicon frame
%> @author Luke Sy (UNSW GSBME)
%> @date 24 Sept 2018
%>
%> @param obj this XsensBody
%> @param qR XsensBody containing transformation quaternion (1 x 4) 
%>           from world frame to vicon frame
%>
%> @retval out XsensBody in vicon frame.
% ======================================================================
function out = toViconFrame(obj, qR)
    out = obj.copy();
    out.frame = 'vicon';
           
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n))
            out.(n).ori = quatmultiply(qR.(n).ori, obj.(n).ori);
        end
    end
end