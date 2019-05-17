function [ vechat ] = hat( vec )
% HAT builds the 3x3 skew symmetric matrix from the 3x1 input or 4x4 from 6x1 input
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implements equations (4) or (5), depending on input size, in the paper.
%
% input: 
%   vec: 3x1 vector phi or 6x1 vector xi
%
% output: 
%   vechat: the 3x3 skew symmetric matrix that can be used
%             to implement the cross product, or 4x4 for transformation
%             matrices
%

validateattributes(vec,{'double'},{'ncols',1})

if size(vec,1) == 3
    
    vechat = [  0,     -vec(3),  vec(2);
            vec(3),   0    , -vec(1);
           -vec(2),  vec(1),   0    ];  
       
elseif size(vec,1) == 6
    vechat = [ hat( vec(4:6,1) ) vec(1:3,1); zeros(1,4) ];    
else   
    warning('hat.m:  invalid vector length for hat operator\n');   
end

