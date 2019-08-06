function [ J ] = vec2jac( vec )
% VEC2JAC Construction of the 3x3 J matrix or 6x6 J matrix
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (98) or (100), depending on the input size, in the paper.
%
% input:
%   vec: a 3x1 vector or 6x1 vector xi
%
% output:
%   J: the 3x3 jacobian matrix or 6x6 jacobian matrix, depending on the input
%

validateattributes(vec,{'double'},{'ncols',1})

tolerance = 1e-12;

if size(vec,1) == 3
    
    phi = vec;
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        J = vec2jacSeries(phi,10);
    else
        axis = phi/ph;

        cph = (1 - cos(ph))/ph;
        sph = sin(ph)/ph;

        J = sph * eye(3) + (1 - sph) * axis * axis' + cph * hat(axis);
    end       
    
elseif size(vec,1) == 6
        
    rho = vec(1:3);
    phi = vec(4:6);
    
    ph = norm(phi);
    if ph < tolerance;
        % If the angle is small, fall back on the series representation
        J = vec2jacSeries(phi,10);
    else
        Jsmall = vec2jac( phi );
        Q = vec2Q( vec );
        J = [ Jsmall Q; zeros(3) Jsmall ];
    end
    
else   
    warning('vec2jac.m:  invalid input vector length\n');   
end

end




