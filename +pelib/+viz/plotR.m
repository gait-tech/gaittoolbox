function plotR(R, origin, color)
	% Plot rotation matrix column basis
	% 
	% :param R: rotation matrix
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 10/11/18

    if nargin <= 1, origin = [0 0 0]; end
    if nargin <= 2, color = 'rgb'; end
    
    for i=1:3
        quiver3(origin(1), origin(2), origin(3), ...
                R(1,i), R(2,i), R(3,i), ...
                'Color', color(i));
    end
end