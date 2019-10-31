function [ T ] = tranValidate( T )
% VALIDATETRANSFORMATION causes an error if the transformation matrix is not in SE(3)
%  
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
% input:
%   T: a 4x4 matrix
%
% output:
%   T: If T is a valid transformation matrix, it is returned.
%      Otherwise, an error is thrown.
%
validateattributes(T,{'double'},{'size',[4,4]});
rotValidate(T(1:3,1:3));

err = max(abs( [0,0,0,1] - T(4,:) ) );
if err > 1e-10
    error('tranValidate.m:  the bottom row of the transformation matrix should be [0,0,0,1]. Error: %f', err)
end

end

