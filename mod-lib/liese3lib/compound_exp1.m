% COMPOUND_EXP1 shows the effect of compounding several uncertain transformations in a row
%
% From: Timothy D Barfoot and Paul T Furgale, 
%       Associating Uncertainty with Three-Dimensional Poses for use in Estimation Problems
%		DOI: 10.1109/TRO.2014.2298059
%       tim.barfoot@utoronto.ca, paul.furgale@mavt.ethz.ch
%
%		Generates figure 2 in the paper.
%
% input: 
%   none - but can change parameters below
%
% output: 
%   compound_exp1.eps : figure showing results
%
clear all;

% Define a PDF over transformations (mean and covariance)
xibar = [1 0 0 0 0 0]';
Tbar = vec2tran( xibar );
cholSigma = diag( [ 0 0 0 0 0 0.03 ] );
Sigma = cholSigma*cholSigma';

% Generate some random samples of the rotation matrices
imax = 1000;
kmax = 101;
for i = 1 : imax
    T{i,1} = vec2tran( zeros(6,1) );
    for k = 2:kmax
        T{i,k} = vec2tran( cholSigma * randn(6,1) )*Tbar * T{i,k-1};
    end
end

% Propagate the uncertainty using second- and fourth-order methods
T_est{1} = vec2tran( zeros(6,1) );
Sigma_est1{1} = zeros(6);  % second order covariance
Sigma_est2{1} = zeros(6); % fourth order covariance
for k = 2:kmax
       
   % Second-order method
   [Ttemp,Sigmatemp] = compound( Tbar, Sigma, T_est{k-1}, Sigma_est1{k-1}, 1 );
   T_est{k} = Ttemp;
   Sigma_est1{k} = Sigmatemp;

   % Fourth-order method
   [Ttemp,Sigmatemp] = compound( Tbar, Sigma, T_est{k-1}, Sigma_est2{k-1}, 2 );
   Sigma_est2{k} = Sigmatemp;

end

% Now plot the transformations
hf=figure(1)
clf
hold on

% Plot some dots that will be overwritten to work on the legend
set(plot( 0, 0),'Color',[0 0 1],'Marker','.','MarkerSize',7,'LineWidth',1);
set(plot( 0, 0),'Color',[1 0 0],'Marker','.','MarkerSize',7,'LineWidth',1);
set(plot( 0, 0),'Color',[0 1 0],'Marker','.','MarkerSize',7,'LineWidth',1);
legend('samples','second-order','fourth-order','Location','SouthWest')

% Plot the random samples' trajectory lines (in a frame attached to the start)
for i = 1 : imax
    for k = 1 : kmax
       r = T{i,k}(1:3,1:3)'*T{i,k}(1:3,4);
       x(k) = r(1);
       y(k) = r(2);
    end
    set(plot( x, y),'Color',[0.8 0.8 1]);
end 

% Plot the random samples' xy-locations (in a frame attached to the start)
kstep = 100;
sigma = 3;
for k = 1 : kstep: kmax
    x = zeros(1,imax);
    y = zeros(1,imax);

    % Plot blue dots for random samples 
    for i = 1 : imax
       r = T{i,k}(1:3,1:3)'*T{i,k}(1:3,4);
       x(i) = r(1);
       y(i) = r(2);
       set(plot( x(i), y(i)),'Color',[0.4 0.4 1],'Marker','.','MarkerSize',7);
    end
    
    % Plot the mean of the samples
    xmean = mean(x);
    ymean = mean(y);
    set(plot( xmean, ymean),'Color',[0 0 1],'Marker','.','MarkerSize',7);
    
    % Plot the covariance of the samples
    vSigma = zeros(2);
    for i = 1:imax
       vSigma = vSigma + [x(i) - xmean; y(i) - ymean ]*[x(i) - xmean; y(i) - ymean ]';
    end
    vSigma = vSigma/imax;
    
    [V,D] = eig(vSigma);
    a = sigma*sqrt( D(1,1) );
    b = sigma*sqrt( D(2,2) );
    mmax = 50;
    clines = zeros(2,mmax); 
    for m = 1 : mmax;   
        v = -pi + 2*pi*(m-1)/(mmax-1);
        p = a*squeeze(V(:,1))*sin(v) + b*squeeze(V(:,2))*cos(v);
        clines(1,m) = p(1) + xmean;
        clines(2,m) = p(2) + ymean;
    end
    set(plot( clines(1,:), clines(2,:) ),'Color',[0.4 0.4 1],'LineWidth',1);
    
    % Plot the propagated mean in SE(3) - projected onto x,y
    r = T_est{k}(1:3,1:3)'*T_est{k}(1:3,4);
    set(plot( r(1), r(2)),'Color',[0 0.6 0],'Marker','.','MarkerSize',7);
    
    % Plot the propagated covariance in SE(3) for second-order method - projected onto x,y
    [V,D] = eig(Sigma_est1{k});
    [Y,I] = sort(diag(D),'descend');
    a = sigma*sqrt( D(I(1),I(1)) );
    b = sigma*sqrt( D(I(2),I(2)) );
    c = sigma*sqrt( D(I(3),I(3)) );
    mmax = 50;
    for n = 1:3
        clines = zeros(2,mmax); 
        for m = 1 : mmax;   
            v = -pi + 2*pi*(m-1)/(mmax-1);
            if n == 1
                p = a*squeeze(V(:,I(1)))*sin(v) + b*squeeze(V(:,I(2)))*cos(v);
            elseif n == 2
                p = b*squeeze(V(:,I(2)))*sin(v) + c*squeeze(V(:,I(3)))*cos(v); 
            elseif n == 3
                p = a*squeeze(V(:,I(1)))*sin(v) + c*squeeze(V(:,I(3)))*cos(v);
            end    
            Ttemp = vec2tran( p )*T_est{k};
            r = Ttemp(1:3,1:3)'*Ttemp(1:3,4);
            clines(1,m) = r(1);
            clines(2,m) = r(2);
        end
        set(plot( clines(1,:), clines(2,:) ),'Color',[1 0 0],'LineWidth',1);
    end
    
    % Plot the propagated covariance in SE(3) for fourth-order method - projected onto x,y
    [V,D] = eig(Sigma_est2{k});
    [Y,I] = sort(diag(D),'descend');
    a = sigma*sqrt( D(I(1),I(1)) );
    b = sigma*sqrt( D(I(2),I(2)) );
    c = sigma*sqrt( D(I(3),I(3)) );
    mmax = 50;
    for n = 1:3
        clines = zeros(2,mmax); 
        for m = 1 : mmax;   
            v = -pi + 2*pi*(m-1)/(mmax-1);
            if n == 1
                p = a*squeeze(V(:,I(1)))*sin(v) + b*squeeze(V(:,I(2)))*cos(v);
            elseif n == 2
                p = b*squeeze(V(:,I(2)))*sin(v) + c*squeeze(V(:,I(3)))*cos(v); 
            elseif n == 3
                p = a*squeeze(V(:,I(1)))*sin(v) + c*squeeze(V(:,I(3)))*cos(v);
            end    
            Ttemp = vec2tran( p )*T_est{k};
            r = Ttemp(1:3,1:3)'*Ttemp(1:3,4);
            clines(1,m) = r(1);
            clines(2,m) = r(2);
        end
        set(plot( clines(1,:), clines(2,:) ),'Color',[0 1 0],'LineWidth',1);
    end
end 

axis equal

set(xlabel('$x$'),'Interpreter','Latex','FontSize',16);
set(ylabel('$y$'),'Interpreter','Latex','FontSize',16);

print -depsc compound_exp1.eps
saveas(hf, 'compound_exp1.fig')


