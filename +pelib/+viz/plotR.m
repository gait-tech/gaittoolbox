% ======================================================================
%> @brief Plot rotation matrix column basis
%> @author Luke Sy
%> @date 11 Oct 2018
%>
%> @param R rotation matrix
% ======================================================================
function plotR(R, origin, color)
    if nargin <= 1, origin = [0 0 0]; end
    if nargin <= 2, color = 'rgb'; end
    
    for i=1:3
        quiver3(origin(1), origin(2), origin(3), ...
                R(1,i), R(2,i), R(3,i), ...
                'Color', color(i));
    end
end