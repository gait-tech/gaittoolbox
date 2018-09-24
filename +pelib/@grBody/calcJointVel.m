% ======================================================================
%> @brief Calculate the joint velocity of grBody
%> @author Luke Sy (UNSW GSBME)
%> @date 24 Sept 2018
%>
%> Example: out = obj.calcJointVel(100, {'LTIO', 'RFEO'})
%> out = struct('LTIO', n x 3, 'RFEO', n x 3)
%>
%> @param obj this grBody
%> @param pts cell array of joints to be calculated (if blank: all joints)
%>
%> @retval out struct of velocities
% ======================================================================
function out = calcJointVel(obj, pts)
    if nargin <= 1
        pts = obj.posList;
    end
    
    out = {}; fs = obj.fs;
    for i=1:length(pts)
        n = pts{i};
        vel = diff(obj.(n), 1, 1)*fs;
        out.(n) = [vel(1,:); vel];
    end
end