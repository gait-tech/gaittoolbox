% ======================================================================
%> @brief set marker points
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> @param subj subject string
%> @param marker marker string
%> @param p n x 3 marker points
%> @param e n x 1 enable boolean vector
% ======================================================================
function p = setPoints(subj, marker, p, e)
    if nargin <=3
        n = size(p, 1);
        e = true(n, 1);
    end
    vicon = ViconNexus;
    vicon.SetTrajectory(subj, marker, p(:,1), p(:,2), p(:,3), e);
end