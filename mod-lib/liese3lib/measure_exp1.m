% MEASURE_EXP1 numerically examines passing uncertain pose and point variables through a nonlinear stereo camera model
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates figures 7 and 8 in the paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   measure_exp1a.eps : figure showing results 
%   measure_exp1b.eps : figure showing results 
%

clear all;

% Define stereo camera matrix
M = 200*[1 0 0  0.25/2; ...
        0 1 0    0; ...
        1 0 0 -0.25/2; ...
        0 1 0    0 ];

% Dilation matrix to perturb landmarks
Q = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% Define pose, landmark, and uncertainty
T = vec2tran( [0;0;0;0;0;0] );
p = [10;10;10;1];

% Compute noisefree measurement
ylin = cam( M, T*p );

% Noise scaling values to loop over
alpha = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0];

mmax = size(alpha,2);
for m=1:mmax

   m

   % Define covariance matrix
   Xi = diag(alpha(m)*[0.1 0.1 0.01 0.01 0.01 0.01 1 1 1]);

   % Compute the expected measurement using Monte Carlo
   ymc = zeros(4,1);
   nsamplemax=1000000;
   sqrtXi = sqrt(Xi);
   for n =1:nsamplemax
      % Draw a sample from the pose/landmark uncertainty
      xi = sqrtXi*randn(9,1);
      Tnoisy = vec2tran( xi(1:6) )*T;
      pnoisy = p + Q*xi(7:9);

      % Compute the measurement
      ymcsample{n} = cam( M, Tnoisy*pnoisy );

      % Compute the mean
      ymc = ymc + ymcsample{n};
   end
   ymc = ymc / nsamplemax;
   ymccov = zeros(4,4);
   for n=1:nsamplemax
      ymccov = ymccov + (ymcsample{n}-ymc)*(ymcsample{n}-ymc)';
   end
   ymccov = ymccov / nsamplemax;

   % Compute the expected measurement using the new technique in the paper
   ynew = ylin;
   I = eye(4);
   F = camJac(M,T*p);
   G = [ point2fs(T*p) T*Q ];
   for j = 1:4
      He{j} = G'*camHess(M,T*p,j)*G;
      for i=1:4
         He{j} = He{j} + F(j,i)*[ point2sf(I(:,i))*point2fs(T*p)  point2sf(I(:,i))*T*Q; ...
                                       (point2sf(I(:,i))*T*Q)'        zeros(3) ];                        
      end           
      ynew = ynew + 0.5*trace(He{j}*Xi)*I(:,j); 
   end
   ynewcov = F*G*Xi*G'*F' - (ynew-ylin)*(ynew-ylin)';
   T4 = zeros(4,4);
   for i=1:4
      for j=1:4
         for kk=1:9
            for ll=1:9
               for mm=1:9
                  for nn=1:9
                     T4(i,j) = T4(i,j) + He{i}(kk,ll)*He{j}(mm,nn) ...
                         * (Xi(kk,ll)*Xi(mm,nn)+Xi(kk,mm)*Xi(ll,nn)+Xi(kk,nn)*Xi(ll,mm));
                  end
               end
            end
         end
      end
   end
   ynewcov = ynewcov + T4/4;

   % Covariance of the linearization approach
   ylincov = F*G*Xi*G'*F';

   % Compute the expected measurement using the sigmapoint transformation
   kappa = 0;
   L = 9;
   S = sqrt(Xi);
   spoint = [zeros(L,1) sqrt(L+kappa)*S -sqrt(L+kappa)*S];
   yspsample{1} = ylin;
   ysp = kappa/(kappa+L)*yspsample{1};
   for n=2:2*L+1
      Tsample = vec2tran( spoint(1:6,n) )*T;
      psample = p + Q*spoint(7:9,n);
      yspsample{n} = cam( M, Tsample*psample );
      ysp = ysp + yspsample{n}/(2*(L+kappa));
   end  
   yspcov = kappa/(kappa+L)*(yspsample{1}-ysp)*(yspsample{1}-ysp)';
   for n=2:2*L+1
      yspcov = yspcov + (yspsample{n}-ysp)*(yspsample{n}-ysp)'/(2*(L+kappa));
   end

   % Plot the left-image means/covariances on the last iteration
   % (largest noise value)
   if m==mmax
       hf=figure(1);
       clf;
       axis equal
       hold on;
       plot( ymc(1), ymc(2), 'b.' );
       plot( ysp(1), ysp(2), 'mo' );
       plot( ylin(1), ylin(2), 'rx' );
       plot( ynew(1), ynew(2), 'gs' );
       set(legend('Monte Carlo','sigmapoint','first/second-order','second/fourth-order', 'Location','NorthWest'),'FontSize',12);

       for n=1:min(nsamplemax,15000)  % don't want plot too crowded
          plot(ymcsample{n}(1),ymcsample{n}(2));
       end
       plot( ymc(1), ymc(2), 'b.' );
       plotcov( ymc(1:2), ymccov(1:2,1:2),1,'b-');

       plot( ysp(1), ysp(2), 'mo' );
       plotcov( ysp(1:2), yspcov(1:2,1:2),1,'m-');

       plot( ylin(1), ylin(2), 'rx' );
       plotcov( ylin(1:2), ylincov(1:2,1:2),1,'r-');

       plot( ynew(1), ynew(2), 'gs' );
       plotcov( ynew(1:2), ynewcov(1:2,1:2),1,'g-');

       axis([100 300 100 300]);

       %set(title('Left Image Feature Uncertainty (one-standard-deviation ellipses)'),'FontSize',14);
       set(xlabel('Horizontal Image Coordinate [pixels]'),'FontSize',14);
       set(ylabel('Vertical Image Coordinate [pixels]'),'FontSize',14);

       print -depsc measure_exp1a.eps
       saveas(hf, 'measure_exp1a.fig')

   end

   elin(m) = sqrt((ylin-ymc)'*(ylin-ymc));
   elincov(m) = sqrt(trace((ylincov-ymccov)'*(ylincov-ymccov)));

   enew(m) = sqrt((ynew-ymc)'*(ynew-ymc));
   enewcov(m) = sqrt(trace((ynewcov-ymccov)'*(ynewcov-ymccov)));

   esp(m) = sqrt((ysp-ymc)'*(ysp-ymc));
   espcov(m) = sqrt(trace((yspcov-ymccov)'*(yspcov-ymccov))); 

end

hf=figure(2);
clf;

subplot(1,2,1);
hold on;
plot( alpha, esp, 'mo-' );
plot( alpha, elin, 'rx-' );
plot( alpha, enew, 'gs-' );
set(xlabel('Noise scaling, \alpha'),'FontSize',14);
set(ylabel('Mean error (compared to Monte Carlo) [pixels], \epsilon_{mean}'),'FontSize',14);
set(legend('sigmapoint','first-order','second-order', 'Location','NorthWest'),'FontSize',12);

subplot(1,2,2);
hold on;
plot( alpha, espcov, 'mo-' );
plot( alpha, elincov, 'rx-' );
plot( alpha, enewcov, 'gs-' );
set(xlabel('Noise scaling, \alpha'),'FontSize',14);
set(ylabel('Covariance error (compared to Monte Carlo) [pixels^2], \epsilon_{cov}'),'FontSize',14);
set(legend('sigmapoint','second-order','fourth-order', 'Location','NorthWest'),'FontSize',12);

print -depsc measure_exp1b.eps
saveas(hf, 'measure_exp1b.fig')
