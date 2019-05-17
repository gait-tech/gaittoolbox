function [ veccurlyhat ] = curlyhat( vec )
% CURLYHAT builds the 6x6 curly hat matrix from the 6x1 input
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implements equation (12) in the paper.  
%
% input: 
%   vec: 6x1 vector xi
%
% output: 
%   vechat: the 6x6 curly hat matrix 
%

validateattributes(vec,{'double'},{'size',[6,1]});

phihat = hat( vec(4:6) );
veccurlyhat = [ phihat hat(vec(1:3)); zeros(3) phihat ];

end
