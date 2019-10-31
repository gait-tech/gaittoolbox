function [ vechat ] = point2fs( vec )
% POINT2FS turns a 4x1 homogeneous point into a special 4x6 matrix
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implements equation (72), left, in the paper.
%
% input: 
%   vec: 4x1 homogeneous point
%
% output: 
%   vechat: the 4x6 matrix used to swap places with pose perturbations
%

validateattributes(vec,{'double'},{'size',[4,1]});

vechat = [ vec(4)*eye(3)  -hat(vec(1:3)); ...
             [0 0 0]         [0 0 0] ];

end