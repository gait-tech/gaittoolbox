function [ Q ] = vec2Q( vec )
% VEC2Q Construction of the 3x3 Q matrix
%   
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Implementation of equation (102), in the paper.
%
% input:
%   vec: a 6x1 vector 
%
% output:
%   Q: the 3x3 Q matrix
%

validateattributes(vec,{'double'},{'size',[6,1]});

rho = vec(1:3);
phi = vec(4:6);

ph = norm(phi);
ph2 = ph*ph;
ph3 = ph2*ph;
ph4 = ph3*ph;
ph5 = ph4*ph;

cph = cos(ph);
sph = sin(ph);

rx = hat(rho);
px = hat(phi);

t1 = 0.5 * rx;
t2 = ((ph - sph)/ph3) * (px*rx + rx*px + px*rx*px);
m3 = ( 1 - 0.5 * ph2 - cph ) / ph4;
t3 = -m3 * (px*px*rx + rx*px*px - 3*px*rx*px);
m4 = 0.5 * ( m3 - 3*(ph - sph - ph3/6)/ph5 );
t4 = -m4 * (px*rx*px*px + px*px*rx*px);

Q = t1 + t2 + t3 + t4;

end