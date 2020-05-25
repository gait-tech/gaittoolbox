% FUSION_EXP1 numerically test the pose fusion technique in the paper
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates figure 5 in the paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   fusion_exp1.eps : figure showing results
%   
clear all;

% Define the groundtruth pose
Ttrue = vec2tran( [1 0 0 0 0 pi/6]' );
invTtrue = inv( Ttrue );

% Covariance matrices for individual estimates
alpha = 5;
kmax = 3;
Sigma{1} = alpha*diag([2 1 1 0.1 0.2 0.1]);
Sigma{2} = alpha*diag([1 3 1 0.1 0.1 0.2]);
Sigma{3} = alpha*diag([1 1 5 0.2 0.1 0.1]);

% precompute inverses and Cholesky factors
for k = 1:kmax
    invSigma{k} = inv( Sigma{k} );
    cholSigma{k} = chol( Sigma{k}, 'lower' );
end

% Loop over many cases and average results
nmax=7;
V = zeros(nmax,1);
ep = zeros(nmax,1);
trcov = zeros(nmax,1);

mmax=1000;
for m=1:mmax

   m

   % Generate some noisy measurements
   T{1} = vec2tran( cholSigma{1}*randn(6,1) )*Ttrue;
   T{2} = vec2tran( cholSigma{2}*randn(6,1) )*Ttrue;
   T{3} = vec2tran( cholSigma{3}*randn(6,1) )*Ttrue;

   % Solve for pose using our algorithm
   for n=1:7
       Test = eye(4);
       for i=1:20      % Gauss-Newton iterations
          LHS = zeros(6);
          RHS = zeros(6,1);
          for k=1:kmax
             xik = tran2vec( Test*inv(T{k}) );
             if n <= 6
                invJ = vec2jacInvSeries( xik, n );
             else
                invJ = vec2jacInv( xik );
             end
             invJtS = invJ'*invSigma{k};
             LHS = LHS + invJtS*invJ;
             RHS = RHS + invJtS*xik;
          end
          xi = -LHS \ RHS;
          Test = vec2tran( xi )*Test;
       end

       % How low did the objective function get?
       for k=1:kmax
          xik = tran2vec( Test*inv(T{k}) );
          V(n) = V(n) + xik'*invSigma{k}*xik/2;
       end

       % How close is the estimate to the true pose?
       xi = tran2vec( Ttrue*inv(Test) );
       ep(n) = ep(n) + xi'*xi;

       % How big is the covariance?
       Sigma_est = inv(LHS);
       trcov(n) = trcov(n) + trace(Sigma_est'*Sigma_est);

   end
end

V = V/mmax
ep = sqrt(ep/mmax)
trcov = sqrt(trcov/mmax)

nterms = [1 2 3 4 5 6 7];

hf=figure(2);
clf;
hold on;
h=subplot(1,2,1);
set(plot( nterms, V, 'ks-'),'LineWidth',1);
set(h, 'XTick', nterms, 'XLim', [1 7] );
set(h,'XTickLabel', ['1' '2' '3' '4' '5' '6' ' ' ]');
ylimits = get(h,'YLim');
text(6.8,ylimits(1)-0.025*(ylimits(2)-ylimits(1)),'$\infty$','interpreter','latex','clipping','off');
set(xlabel('Number of terms, N'),'FontSize',14);
set(ylabel('Average final cost function, V'),'FontSize',14);

h=subplot(1,2,2);
set(plot( nterms, ep, 'ks-' ),'LineWidth',1);
set(h, 'XTick', nterms, 'XLim', [1 7] );
set(h,'XTickLabel', ['1' '2' '3' '4' '5' '6' ' ' ]');
ylimits = get(h,'YLim');
text(6.8,ylimits(1)-0.025*(ylimits(2)-ylimits(1)),'$\infty$','interpreter','latex','clipping','off');
set(xlabel('Number of terms, N'),'FontSize',14);
set(ylabel('RMS error (compared to true pose), \epsilon'),'FontSize',14);

print -depsc fusion_exp1.eps
saveas(hf, 'fusion_exp1.fig')