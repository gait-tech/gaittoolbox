function Hi = camHess( M, g, i )
% CAMHESS produces the Hessian of a stereo camera model
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		See section V.E in the paper (last set of equations).      
%
% input: 
%   M: 4x4 stereo camera matrix
%   g: 4x1 homogeneous point
%	i: scalar to indicate which component of Hessian
%
% output: 
%   Hi: 4x4 Hessian of ith row of stereo camera model
%

validateattributes(M,{'double'},{'size',[4,4]});
validateattributes(g,{'double'},{'size',[4,1]});
validateattributes(i,{'numeric'},{'scalar'});

b = 2*M(1,4)/M(1,1);
if i==1
  Hi = [ 0 0 -1 0; ...
         0 0  0 0; ...
        -1 0 (2*g(1)+b*g(4))/g(3) -b/2; ...
         0 0 -b/2 0]*M(1,1)/g(3)^2; 
elseif i==2
  Hi = [ 0  0  0 0; ...
         0  0 -1 0; ...
         0 -1 2*g(2)/g(3) 0; ...
         0  0  0 0]*M(2,2)/g(3)^2;
elseif i==3
  Hi = [ 0 0 -1 0; ...
         0 0  0 0; ...
        -1 0 (2*g(1)-b*g(4))/g(3) b/2; ...
         0 0 b/2 0]*M(1,1)/g(3)^2;        
elseif i==4
  Hi = [ 0  0  0 0; ...
         0  0 -1 0; ...
         0 -1 2*g(2)/g(3) 0; ...
         0  0  0 0]*M(2,2)/g(3)^2;     
else
  warning('camHess.m: invalid measurement row\n');
end  

end