% ======================================================================
%> @brief get marker points
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> @param subj subject string
%> @param marker marker string
%> @param idx [Optional] index of points to be returned. Default= 1:end
%>
%> @retval p n x 3 marker points
% ======================================================================
function p = getPoints(subj, marker, idx)
    vicon = ViconNexus;
    [x, y, z, e] = vicon.GetTrajectory(subj, marker);
    
    if nargin <= 2
        idx = 1:length(x);
    end
    
    p = [x(idx)' y(idx)' z(idx)'];
end