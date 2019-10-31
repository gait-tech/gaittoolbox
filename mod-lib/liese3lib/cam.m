function h = cam( M, g ) 
% CAM produces stereo pixel measurements given a point in camera coordinates
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		See section V.E in the paper (first equation).      
%
% input: 
%   M: 4x4 stereo camera matrix
%   g: 4x1 homogeneous point
%
% output: 
%   h: 4x1 pixel coords (horiz left, vert left, horiz right, vert right)
%

validateattributes(M,{'double'},{'size',[4,4]});
validateattributes(g,{'double'},{'size',[4,1]});

h = M*g/g(3);

end