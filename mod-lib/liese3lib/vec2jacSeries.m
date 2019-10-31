function [ J ] = vec2jacSeries( vec, N )
% VEC2JACSERIES Construction of the J matrix from Taylor series
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (98) or (100), depending on the input size, in the paper.
%
% input:
%   phi: a 3x1 vector
%   N:   number of terms to include in the series
%
% output:
%   J: the 3x3 J matrix 
%

validateattributes(vec,{'double'},{'ncols',1})
validateattributes(N,{'numeric'},{'scalar'});

if size(vec,1) == 3 
    
    J = eye(3);
    pxn = eye(3);
    px = hat(vec);
    for n = 1:N
        pxn = pxn*px/(n + 1);    
        J = J + pxn;
    end
    
elseif size(vec,1) == 6    
    
    J = eye(6);
    pxn = eye(6);
    px = curlyhat(vec);
    for n = 1:N
        pxn = pxn*px /(n + 1);    
        J = J + pxn;
    end    
    
else   
    warning('vec2jacSeries.m:  invalid input vector length\n');   
end

end

