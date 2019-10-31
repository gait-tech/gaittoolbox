function [ invJ ] = vec2jacInvSeries( vec, N )
% VEC2JACINVSERIES Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix. Series representation
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (99) or (103), depending on the input size, in the paper.
%
% input: 
%   vec: 3x1 vector or 6x1 vector
%   N:   number of terms to include in the series
%
% output: 
%   invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix
%

validateattributes(vec,{'double'},{'ncols',1})
validateattributes(N,{'numeric'},{'scalar'});

if size(vec,1) == 3 
    
    invJ = eye(3);
    pxn = eye(3);
    px = hat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    
elseif size(vec,1) == 6    
    
    invJ = eye(6);
    pxn = eye(6);
    px = curlyhat(vec);
    for n = 1:N
        pxn = pxn * px/n;
        invJ = invJ + bernoullinumber(n) * pxn;
    end
    
else   
    warning('vec2jacInvSeries.m:  invalid input vector length\n');   
end

end

