% MEASURE_EXP2 shows the effect of an uncertain pose applied to a constant vector in 3D
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates a bonus figure not shown in the final paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   measure_exp2.eps : figure showing results 
%

% Define a PDF over rotations (mean and covariance)
phibar = [0 0 pi/6]';
Cbar = vec2rot( phibar );
R = vec2rot( [pi/6 pi/4 0]' );
Sigma = Cbar*R*diag( [ 0.01^2 0.1^2 0.6^2] )*R'*Cbar';
cholSigma = chol( Sigma,'lower' );

% Generate some random samples of the rotation matrices
imax = 1000;
for i = 1 : imax    
    C{i} = vec2rot( cholSigma * randn(3,1) )*Cbar;
end

% Now try to compute the mean
Cbar_ran = eye(3);
for j = 1 : 10
    del = zeros(3,1);
    for i = 1 : imax
        del = del + rot2vec( C{i}*Cbar_ran' );
    end
    del = del / imax;
    Cbar_ran = vec2rot( del ) * Cbar_ran;
    rotValidate(Cbar_ran);
end
Cbar_ran

% Generate some sigmapoint samples of the rotation matrix
Csp{1} = Cbar;
pmax = 3;
kappa = 0;
alpha = sqrt(pmax+kappa);
for p = 1 : pmax
   Csp{1+p} = vec2rot( alpha * squeeze(cholSigma(:,p) ))*Cbar;
   Csp{1+p+pmax} = vec2rot( -alpha * squeeze(cholSigma(:,p) ))*Cbar;
end

% Now try to compute the mean
Cbar_sp = eye(3);
for j = 1 : 10
    del = kappa/(pmax+kappa)*rot2vec( Csp{1}*Cbar_sp' );
    for p = 2 : 2*pmax+1
        del = del + 1/(2*(pmax+kappa))*rot2vec( Csp{p}*Cbar_sp' );
    end
    Cbar_sp = vec2rot( del ) * Cbar_sp;
    rotValidate(Cbar_sp);
end
Cbar_sp

% Now plot the rotation times a constant vector
figure(1)
clf
hold on

% To visualize, plot a rotated vector several different ways:
x = [1 0 0]';

% Vector sample mean (average of sampled rotations times vector x)
ybar = zeros(3,1);
for i = 1 : imax
    ybar = ybar + C{i} * x;
end
ybar = ybar / imax
set(plot3([0 ybar(1)], [0 ybar(2)], [0 ybar(3)] ),'Color','g', 'LineWidth',2,'Marker','.','MarkerSize', 7);

% Vector mean via sigmapoint method
ysp = kappa/(pmax+kappa) * Csp{1} * x;
for p = 2: 2*pmax+1
    ysp = ysp + 1/(2*(pmax+kappa)) * Csp{p} * x;
end
set(plot3([0 ysp(1)], [0 ysp(2)], [0 ysp(3)] ),'Color','m', 'LineWidth',2,'Marker','.','MarkerSize', 7);

% Vector mean via second-order Taylor-series expansion (should also be too short)
ytay2 = (eye(3) - trace(Sigma)*eye(3)/2 + Sigma/2 )*Cbar*x
set(plot3([0 ytay2(1)], [0 ytay2(2)], [0 ytay2(3)] ),'Color',[0.91 0.41 0.17], 'LineWidth',2,'Marker','.','MarkerSize', 7);

% Vector mean via fourth-order Taylor-series expansion (should also be too short)
ytay4 = ( eye(3) - trace(Sigma)*eye(3)/2 + Sigma/2 + ((trace(Sigma))^2 + 2*trace(Sigma^2))*eye(3)/24 - Sigma*(trace(Sigma)*eye(3) + 2*Sigma)/24 )*Cbar*x
set(plot3([0 ytay4(1)], [0 ytay4(2)], [0 ytay4(3)] ),'Color','c', 'LineWidth',2,'Marker','.','MarkerSize', 7);


% Groundtruth (Cbar times vector x)
ytrue = Cbar*x
set(plot3([0 ytrue(1)], [0 ytrue(2)], [0 ytrue(3)] ),'Color','k','LineWidth',2,'Marker','.','MarkerSize', 7);

% Legend/
legend('sampled','sigmapoint','2nd-order','4th-order','noisefree','Location','NorthWest');
set(xlabel('$x$'),'Interpreter','Latex','FontSize',16);
set(ylabel('$y$'), 'Interpreter','Latex','FontSize',16);
set(zlabel('$z$'), 'Interpreter','Latex','FontSize',16);

% Plot equiprobable contours - three great circles of original uncertainty
% ellipsoid mapped to the vector space
kmax = 3;
nmax = 3;
mmax = 51;
[V,D] = eig(Sigma);
a = sqrt( D(1,1) );
b = sqrt( D(2,2) );
c = sqrt( D(3,3) );
clines = zeros(3,3,nmax,mmax); 
for n = 1 : nmax
    for m = 1 : mmax;   
        v = -pi + 2*pi*(m-1)/(mmax-1);
        if n == 1
            delphi = b*squeeze(V(:,2))*sin(v) + c*squeeze(V(:,3))*cos(v);
        elseif n == 2
            delphi = a*squeeze(V(:,1))*sin(v) + c*squeeze(V(:,3))*cos(v);
        elseif n == 3
            delphi = a*squeeze(V(:,1))*sin(v) + b*squeeze(V(:,2))*cos(v);
        end
        for k = 1 : kmax  
            y = vec2rot( k*delphi )*Cbar*x;
            clines(k,1,n,m) = y(1);
            clines(k,2,n,m) = y(2);
            clines(k,3,n,m) = y(3);
        end
    end
end
for n=1:nmax
    for k=1:kmax
        set(plot3( squeeze(clines(k,1,n,:)), squeeze(clines(k,2,n,:)), squeeze(clines(k,3,n,:)) ),'Color',[(k-1)/2.5 (k-1)/2.5 1],'LineWidth',1);
    end
end

% Plot the random samples
for i = 1 : min(1000,imax)
    y = C{i} * x;
    set(plot3( y(1), y(2), y(3) ),'Color',[0.4 0.4 1],'Marker','.','MarkerSize', 7);
end

% Plot the sigmapoint samples
for p = 1:2*pmax+1
    y = Csp{p} * x;
    set(plot3( y(1), y(2), y(3) ),'Color','m','Marker','.','MarkerSize', 10);
end   

axis equal
axis([-0.1 1 -0.8 0.85 -0.6 1])

print -depsc measure_exp2.eps
saveas(hf, 'measure_exp2.fig')

