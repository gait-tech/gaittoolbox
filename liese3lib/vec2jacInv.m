function [ invJ ] = vec2jacInv( vec )
% VEC2JACINV Construction of the 3x3 J^-1 matrix or 6x6 J^-1 matrix in closed form
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
%
% output: 
%   invJ: 3x3 inv(J) matrix or 6x6 inv(J) matrix
%
   
validateattributes(vec,{'double'},{'ncols',1})

tolerance = 1e-12;

if size(vec,1) == 3
    
    phi = vec;
    
    ph = norm(phi);
    if ph < tolerance
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(phi,10);
    else
        axis = phi/norm(phi);
        ph_2 = 0.5*ph;

        invJ =   ph_2 * cot(ph_2)* eye(3)...
               + (1 - ph_2 * cot(ph_2))* axis * axis'...
               - ph_2 * hat(axis);
    end   
    
elseif size(vec,1) == 6

    rho = vec(1:3);
    phi = vec(4:6);
    
    ph = norm(phi);
    if ph < tolerance;
        % If the angle is small, fall back on the series representation
        invJ = vec2jacInvSeries(phi,10);
    else
        invJsmall = vec2jacInv( phi );
        Q = vec2Q( vec );
        invJ = [ invJsmall -invJsmall*Q*invJsmall; zeros(3) invJsmall ];
    end
    
else   
    warning('vec2jacInv.m:  invalid input vector length\n');   
end   
   
end

