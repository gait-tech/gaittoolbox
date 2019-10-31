function g = camInv( M, h ) 
% CAMINV produces a point in camera coordinates from stereo pixel measurements
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		See section V.E in the paper (inverse of first equation).      
%
% input: 
%   M: 4x4 stereo camera matrix
%   h: 4x1 pixel coords (horiz left, vert left, horiz right, vert right)
%
% output: 
%   g: 4x1 homogeneous point
%

validateattributes(M,{'double'},{'size',[4,4]});
validateattributes(h,{'double'},{'size',[4,1]});

g = zeros(4,1);

g(4) = 1.0;
g(3) = 2*M(1,4)/(h(1)-h(3))*g(4);
g(1) = g(3)*(h(1)+h(3)-2*M(1,3))/(2*M(1,1));
g(2) = g(3)*(h(2)-M(2,3))/M(2,2);

end