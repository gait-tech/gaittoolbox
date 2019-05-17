function H = camJac( M, g ) 
% CAMJAC produces the Jacbian of a stereo camera model
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		See section V.E in the paper (third displayed equation).      
%
% input: 
%   M: 4x4 stereo camera matrix
%   g: 4x1 homogeneous point
%
% output: 
%   H: 4x4 Jacobian
%

validateattributes(M,{'double'},{'size',[4,4]});
validateattributes(g,{'double'},{'size',[4,1]});

g3i = 1/g(3);
g3is = g3i*g3i;
H = M*[ g3i   0  -g(1)*g3is   0; ... 
         0   g3i -g(2)*g3is   0; ...
         0    0       0       0; ...
         0    0  -g(4)*g3is  g3i];
    
end