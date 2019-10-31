% FUSION_EXP2 numerically test the pose fusion technique in the paper
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates figure 6 in the paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   fusion_exp2.eps : figure showing results
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

% precompute inverses, Cholesky factors, and noisy measurements
for k = 1:kmax
    invSigma{k} = inv( Sigma{k} );
    cholSigma{k} = chol( Sigma{k}, 'lower' );
    T{k} = vec2tran( cholSigma{k}*randn(6,1) )*Ttrue;
end

% Loop over many cases and average results
imax=15;
nmax=3;
V = zeros(imax,nmax);

for n=1:3   % different order of algorithms
    Test = eye(4);
    for i=1:imax      % Gauss-Newton iterations

     % How low did the objective function get?
      V(i,n) = 0;
      for k=1:kmax
        xik = tran2vec( Test*inv(T{k}) );
        V(i,n) = V(i,n) + xik'*invSigma{k}*xik/2;
      end

      % Solve for the optimal pose change and apply it
      LHS = zeros(6);
      RHS = zeros(6,1);
      for k=1:kmax
         xik = tran2vec( Test*inv(T{k}) );
         if n==1
            invJ = vec2jacInvSeries( xik, 2 );  % second order
         elseif n==2
            invJ = vec2jacInvSeries( xik, 4 );  % fourth order
         elseif n==3
            invJ = vec2jacInv( xik );  % infinite order
         end
         invJtS = invJ'*invSigma{k};
         LHS = LHS + invJtS*invJ;
         RHS = RHS + invJtS*xik;
      end
      xi = -LHS \ RHS;
      Test = vec2tran( xi )*Test;

    end
end

hf=figure(2);
clf;

h=subplot(1,2,1);
hold on;
set(plot( V(:,1), 'rx-'),'LineWidth',1);
set(plot( V(:,2), 'gs-'),'LineWidth',1);
set(plot( V(:,3), 'k.-'),'LineWidth',1);
set(h, 'XTick', linspace(1,imax,imax), 'XLim', [1 imax] );
set(legend('second-order','fourth-order','infinite-order', 'Location','NorthEast'),'FontSize',12);
set(xlabel('Number of Iterations'),'FontSize',14);
set(ylabel('Average final cost function, V'),'FontSize',14);

h=subplot(1,2,2);
hold on;
set(plot( V(:,1), 'rx-'),'LineWidth',1);
set(plot( V(:,2), 'gs-'),'LineWidth',1);
set(plot( V(:,3), 'k.-'),'LineWidth',1);
set(h, 'XTick', linspace(4,imax,imax-4+1), 'XLim', [4 imax] );
a = min(V(imax,:));
b = max(V(imax,:));
d = 1/3*(b-a);
set(h,'YLim',[a-d b+2*d] );
set(xlabel('Number of Iterations'),'FontSize',14);
set(ylabel('Average final cost function, V'),'FontSize',14);

print -depsc fusion_exp2.eps
saveas(hf, 'fusion_exp2.fig')