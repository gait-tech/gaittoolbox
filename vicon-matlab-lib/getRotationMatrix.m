% ======================================================================
%> @brief get marker points
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> Origin is assume at point 1.
%> basis vector X generated from p2 - p1;
%> basis vector Y generated from p3 - p1 - projection with basis vector X
%> basis vector Z generated from cross product of basis vector X and Y
%>
%> @param p1 point 1 (1 x 3 vector assumed as origin)
%> @param p2 point 2 (1 x 3 vector)
%> @param p3 point 3 (1 x 3 vector)
%>
%> @retval R 3 x 3 rotation matrix [basisX basisY basisZ]
% ======================================================================
function R = getRotationMatrix(p1, p2, p3)
    basisX = p2 - p1; 
    basisX = basisX / norm(basisX);
    
    basisY = p3 - p1;
    basisY = basisY - dot(basisY, basisX) * basisX;
    basisY = basisY / norm(basisY);
    
    basisZ = cross(basisX, basisY);
    
    R = [basisX' basisY' basisZ'];
end