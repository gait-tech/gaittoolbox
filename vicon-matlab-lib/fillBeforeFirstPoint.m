% ======================================================================
%> @brief Fill the gap before first point by back fill
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> @param subj subject string
%> @param marker marker string
% ======================================================================
function fillBeforeFirstPoint(subj, marker)
    vicon = ViconNexus;
    [x, y, z, e] = vicon.GetTrajectory(subj, marker);

    firstIdx = max(find(e, 1), 2);
    x(1:firstIdx-1) = x(firstIdx);
    y(1:firstIdx-1) = y(firstIdx);
    z(1:firstIdx-1) = z(firstIdx);
    e(1:firstIdx-1) = e(firstIdx);

    vicon.SetTrajectory(subj, marker, x, y, z, e);
end