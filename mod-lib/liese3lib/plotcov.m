function plotcov( mu, C, k, s )
% PLOTCOV plots a covariance matrix
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
% input: 
%   mu: mean
%   C: covariance matrix
%   k: number of standard deviations
%   s: line type as a string
%
% output: 
%   none
%

   V = zeros(2,2);
   D = zeros(2,2);
   [V, D] = eig(C);
   [a, i] = max([D(1,1), D(2,2)]);
   
   if a == 0
     return
   end
   
   b = min([D(1,1), D(2,2)]);
   
   a = k*sqrt(a);
   b = k*sqrt(b);
   
   beta = atan2( V(2,i), V(1,i) );
   R = [cos(beta) -sin(beta); sin(beta) cos(beta) ];
   
   e = sqrt( 1 - (b/a)^2 );
      
   imax = 1000;
   x = zeros(imax+1,1);
   y = zeros(imax+1,1);
   for i = 1:imax
       theta = 2*pi*(i-1)/(imax-1);
       r = a*(1-e^2)/(1+e*cos(theta)); 
       
       z = R*[r*cos(theta) - (1-e)*a + a; r*sin(theta)] + mu;
       
       x(i) = z(1);
       y(i) = z(2);
   end
   x(imax+1) = x(1);
   y(imax+1) = y(1);
   
   plot( x, y, s);
   
end
