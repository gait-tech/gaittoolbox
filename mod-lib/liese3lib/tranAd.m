function [ AdT ] = tranAd( T )
% TRANAD builds the 6x6 transformation matrix from the 4x4 one
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (11), in the paper.
%
% input: 
%   T: 4x4 transformation matrix
%
% output: 
%   AdT: 6x6 transformation matrix
%

tranValidate( T );

C = T(1:3,1:3);
r = T(1:3,4);
AdT = [ C hat(r)*C; zeros(3) C];

