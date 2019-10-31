function [ T, Sigma ] = compound( T1, Sigma1, T2, Sigma2, method )
% COMPOUND compounds two uncertain transformations
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Primarily implements equations (38) and (55) in the paper (for method=1 or method=2).  
%		See also equation (57) when using method=3 and section III.D for a summary of all
%		four methods.
%
% input: 
%   T1:     4x4 mean of left transformation 
%   Sigma1: 6x6 covariance of left transformation
%   T2:     4x4 mean of right transformation
%   Sigma2: 6x6 covariance of right transformations
%   method: integer indicating method to be used to perform compounding
%           (1=second-order, 2=fourth-order, 3=sigmapoint, 4=Monte Carlo)
%
% output: 
%   T:      4x4 mean of compounded transformation 
%   Sigma:  6x6 covariance of compounded transformation
%

% check the input tranformations are valid
tranValidate( T1 );
tranValidate( T2 );
validateattributes(Sigma1,{'double'},{'size',[6,6]});
validateattributes(Sigma2,{'double'},{'size',[6,6]});

% compound the means
T = T1*T2;
tranValidate( T );

% compute Adjoint of left transformation mean
AdT1 = tranAd( T1 );
Sigma2prime = AdT1 * Sigma2 * AdT1'; 

% compound the covariances based on the chosen method
if method == 1
    
   % Second-order method
   Sigma = Sigma1 + Sigma2prime;

elseif method == 2
   
   % Fourth-order method
   Sigma1rr = Sigma1(1:3,1:3);
   Sigma1rp = Sigma1(1:3,4:6);
   Sigma1pp = Sigma1(4:6,4:6);
   
   Sigma2rr = Sigma2prime(1:3,1:3);
   Sigma2rp = Sigma2prime(1:3,4:6);
   Sigma2pp = Sigma2prime(4:6,4:6);
   
   A1 = [covop1(Sigma1pp) covop1(Sigma1rp+Sigma1rp'); zeros(3) covop1(Sigma1pp) ];
   A2 = [covop1(Sigma2pp) covop1(Sigma2rp+Sigma2rp'); zeros(3) covop1(Sigma2pp) ];
   
   Brr = covop2(Sigma1pp,Sigma2rr)+covop2(Sigma1rp',Sigma2rp)+covop2(Sigma1rp,Sigma2rp')+covop2(Sigma1rr,Sigma2pp);
   Brp = covop2(Sigma1pp,Sigma2rp')+covop2(Sigma1rp',Sigma2pp);
   Bpp = covop2(Sigma1pp,Sigma2pp);
   B = [Brr Brp; Brp' Bpp];
   
   Sigma = Sigma1 + Sigma2prime + (A1*Sigma2prime + Sigma2prime*A1' + A2*Sigma1 + Sigma1*A2')/12 + B/4;

elseif method == 3
        
    % Note 1:  we could get a  speedup in implementation here by noting that 
    % this method is algebraically equivalent to the second-order scheme 
    % (method 1) for the case that the two transformation are uncorrelated
    %
    % Note 2:  if any of the sigmapoints for the rotational components end
    % up being larger than \pi, this method can fail...it's wise to scale
    % down the size of the sigmapoints significantly to avoid this
    
    % Sigmapoint method
    kappa = 1;
    L = 12;
    S = sqrt(kappa)*[ chol( Sigma1, 'lower' ) zeros(6); zeros(6) AdT1*chol( Sigma2, 'lower' ) ];
    
    Sigma = zeros(6);
    for i=1:L
        % pull out sigmapoint for individual transformations
        xi1 = S(1:6,i);
        xi2 = S(7:12,i);
        
        % catch the error concerning making sigmapoints too big
        if norm(xi1(4:6)) > pi || norm(xi2(4:6)) > pi
           error('compound.m:  sigmapoint too large on rotational degree of freedom...try reducing kappa'); 
        end
        
        % positive weight        
        xi = tran2vec( vec2tran(xi1)*vec2tran(xi2) );
        Sigma = Sigma + xi*xi';
        
        % negative weight
        xi = tran2vec( vec2tran(-xi1)*vec2tran(-xi2) );
        Sigma = Sigma + xi*xi';
    end
    Sigma = Sigma / (2*kappa);
    
elseif method == 4
    
    % Monte Carlo method
    nsamples = 1000000;
    invT1 = inv( T1 );
    cholSigma1 = chol( Sigma1, 'lower' );
    cholSigma2 = AdT1*chol( Sigma2, 'lower' );
    
    Sigma = zeros(6);
    for i=1:nsamples
       xi1 = cholSigma1 * randn(6,1);
       xi2 = cholSigma2 * randn(6,1);
       xi = tran2vec( vec2tran(xi1)*vec2tran(xi2) );
       Sigma = Sigma + xi*xi';
    end
    Sigma = Sigma / nsamples;
   
else
   error('compound.m:  not a validate compounding method');
end

end

% covariance operation 1
function A = covop1( B )

validateattributes( B, {'double'},{'size',[3,3]} );

A = -trace(B)*eye(3) + B;

end

% covariance operation 2
function A = covop2( B, C )

validateattributes( B, {'double'},{'size',[3,3]} );
validateattributes( C, {'double'},{'size',[3,3]} );

A = covop1(B)*covop1(C) + covop1(C*B);

end




