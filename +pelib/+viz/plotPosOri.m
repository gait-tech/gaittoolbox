% ======================================================================
%> @brief Plot point with axis
%>
%> @param pos position of point
%> @param ori orientation of segment
%>
% ======================================================================
function plotPosOri(pos, ori)
    hold on;
    for i=1:3
        quiver3(pos(1), pos(2), pos(3), ori(1,i), ori(2,i), ori(3,i));
    end
end