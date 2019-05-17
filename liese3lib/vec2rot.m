function [ C ] = vec2rot( phi )
% VEC2ROT Build a rotation matrix using the exponential map
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (97), in the paper.
%
% input: 
%   phi: 3x1 vector 
%
% output: 
%   C: 3x3 rotation matrix
%

validateattributes(phi,{'double'},{'size',[3,1]});

tolerance = 1e-12;

% Check for a small angle.
% 
angle = norm(phi);
if angle < tolerance
    % If the angle is small, fall back on the series representation.
    % In my experience this is very accurate for small phi
    C = vec2rotSeries(phi,10);
    %display('vec2rot.m:  used series method');
else
    axis = phi/angle;

    cp = cos(angle);
    sp = sin(angle);

    C = cp * eye(3) + (1 - cp) * axis * axis' + sp * hat(axis);
    %display('vec2rot.m:  used analytical method');
end


end

