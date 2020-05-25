function [ T ] = vec2tranSeries( p, N )
% EXPT Build a transformation matrix using the exponential map series with N elements in the series
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (96), series version, in the paper.
%
% input: 
%   p: 6x1 vector
%   N:   number of terms to include in the series
%
% output: 
%   T: 4x4 transformation matrix
% 
%

validateattributes(p,{'double'},{'size',[6,1]});
validateattributes(N,{'numeric'},{'scalar'});

T = eye(4);
xM = eye(4);
bpP = hat(p);
for n = 1:N
    xM = xM * (bpP / n);    
    T = T + xM;
end

end
