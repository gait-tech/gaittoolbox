function [ p ] = tran2vec( T )
% LOGT Compute the matrix log of the transformation matrix T.
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of the algorithm 2 in Appendix B, in the paper.
%
% input:
%   T: a 4x4 transformation matrix
%
% output:
%   p: a 6x1 vector in tangent coordinates computed from T
%

tranValidate(T);

C = T(1:3,1:3);
r = T(1:3,4);

phi = rot2vec(C);

invJ = vec2jacInv(phi);
%invJ = vec2jacInvSeries(phi,10);

rho = invJ * r;
p = [rho;phi];

end

