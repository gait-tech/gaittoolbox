% ======================================================================
%> @brief Fill the gap after last point by forward fill
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> @param subj subject string
%> @param marker marker string
% ======================================================================
function fillAfterLastPoint(subj, marker)
    vicon = ViconNexus;
    [x, y, z, e] = vicon.GetTrajectory(subj, marker);

    lastIdx = find(e, 1, 'last');
    x(lastIdx+1:end) = x(lastIdx);
    y(lastIdx+1:end) = y(lastIdx);
    z(lastIdx+1:end) = z(lastIdx);
    e(lastIdx+1:end) = e(lastIdx);

    vicon.SetTrajectory(subj, marker, x, y, z, e);
end