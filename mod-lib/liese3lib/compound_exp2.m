% COMPOUND_EXP2 numerically compares different methods of compounding two uncertain transformations
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates figure 3 in the paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   compound_exp2.eps : figure showing results
%
clear all;

% Loop over different noise values
imax = 21;
for i=1:imax   

   i

   % Pick two transformations with uncertainties to compound
   if i==1
      alpha(i) = 0.001;
   else
      alpha(i) = (i-1)/(imax-1);
   end
   
   T1 = vec2tran( [0;2;0;pi/6;0;0] );
   Sigma1 = alpha(i)*5*diag([2 1 1 0.1 0.2 0.1]);

   T2 = vec2tran( [0;0;1;pi/4;0;0] );
   Sigma2 = alpha(i)*5*diag([1 2 1 0.1 0.1 0.2]);

   % Work out the answer using our first-order formula
   [Ttemp,Sigma_so] = compound( T1, Sigma1, T2, Sigma2, 1 );

   % Work out the answer using our fourth-order formula
   [Ttemp,Sigma_fo] = compound( T1, Sigma1, T2, Sigma2, 2 );
   
   % Work out the answer using the Sigmapoint Transformation
   [Ttemp,Sigma_sp] = compound( T1, Sigma1, T2, Sigma2, 3 );

   % Work out the answer using Monte Carlo with a large number of samples
   [Ttemp,Sigma_mc] = compound( T1, Sigma1, T2, Sigma2, 4 );
   
   % statistics
   err_so(i) = sqrt( trace( ( Sigma_so - Sigma_mc )'*( Sigma_so - Sigma_mc ) ) );
   err_fo(i) = sqrt( trace( ( Sigma_fo - Sigma_mc )'*( Sigma_fo - Sigma_mc ) ) );
   err_sp(i) = sqrt( trace( ( Sigma_sp - Sigma_mc )'*( Sigma_sp - Sigma_mc ) ) );

end

hf=figure(1);
clf;
hold on;
plot( alpha, err_sp, 'mo-' );
plot( alpha, err_so, 'rx-' );
plot( alpha, err_fo, 'gs-' );
set(xlabel('Noise scaling, \alpha'),'FontSize',16);
set(ylabel('Covariance error (compared to Monte Carlo), \epsilon'),'FontSize',16);
set(legend('sigmapoint','second-order','fourth-order', 'Location','NorthWest'),'FontSize',12);

print -depsc compound_exp2.eps
saveas(hf, 'compound_exp2.fig')


