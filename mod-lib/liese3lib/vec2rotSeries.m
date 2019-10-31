function [ C ] = vec2rotSeries( phi, N )
% VEC2ROTSERIES Build a rotation matrix using the exponential map series with N elements in the series
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (97), in the paper.
%
% input: 3x1 vector phi
%        scalar N, the number of terms to include
%
% output: 3x3 rotation matrix C
%
% input: 
%   phi: 3x1 vector 
%   N:   number of terms to include in the series
%
% output: 
%   C: 3x3 rotation matrix
%

validateattributes(phi,{'double'},{'size',[3,1]});
validateattributes(N,{'numeric'},{'scalar'});

C = eye(3);
xM = eye(3);
cmPhi = hat(phi);
for n = 1:N
    xM = xM * (cmPhi / n);
    C = C + xM;
end

% Project the resulting rotation matrix back onto SO(3)
C = C *inv( sqrtm( C'*C ) );
rotValidate(C);

end

