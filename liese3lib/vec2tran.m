function [ T ] = vec2tran( p )
% EXPT Build a transformation matrix using the exponential map, closed form
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of algorithm 1 in Appendix B, in the paper.
%
% input: 
%   p: 6x1 vector
%
% output: 
%   T: 4x4 transformation matrix
% 

validateattributes(p,{'double'},{'size',[6,1]});

rho = p(1:3);
phi = p(4:6);

C = vec2rot(phi);
J = vec2jac(phi);

T = eye(4);
T(1:3,1:3) = C;
T(1:3,4) = J*rho;

end

