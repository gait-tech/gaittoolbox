function UKF_Rel_LPVA_plot(N_MP,fs,actBody,estBody,pIntS)
%close all
plotWidth = 1000;
plotHeight = 800;
nDim = 3;
% plotMP = 0;
% plotLA = 3;
% plotRA = 6;
N_MP = N_MP-1;
%actState = [actBody.MIDPEL actBody.LTIO actBody.RTIO];
estBodyInt = nan(N_MP,3);
actBodyInt = nan(N_MP,3);

if(strcmp(pIntS,'MP'))
    %pInt = plotMP;
    error('cant plot MP rel to MP')
elseif(strcmp(pIntS,'LA'))
    estBodyInt = estBody.LTIO;
    actBodyInt = actBody.LTIO;
    %pInt = plotLA;
elseif(strcmp(pIntS,'RA'))
    estBodyInt = estBody.RTIO;
    actBodyInt = actBody.RTIO;
    %pInt = plotRA;
else
    error('invalid body part')
end
%% accel and vel plots
%pInt = plotRA;
%Accel Plot
% figure
% for k=1:3                              % plot results
%   subplot(nDim,1,k)
%   hold on
%   vAcc = diff(actState(1:N_MP,k+3+2*pInt))*fs;
%   plot((1:N_MP)/fs, x_rec(k+6+3*pInt,1:N_MP), '-', (1:N_MP)/fs, z_rec(k+pInt,1:N_MP), '--')
%   plot((2:N_MP-1)/fs, vAcc(2:N_MP-1), '.')
%   hold off
%   ylabel(strcat('state ',num2str(k)))
%   xlabel('time (s)')
%   legend('est acc','meas acc', 'tru acc','location','northwest')
% end
% suptitle('accel from imu')


% figure
% hold on
% for k=1:3                              % plot results
%   %subplot(nStates,1,k)
%   plot((1:N_MP-1)/fs, diff(actState(:,k+3))*fs, '.')
% end
% ylabel(strcat('state ',num2str(k)))
% xlabel('time (s)')
% legend('rec pos','tru pos')
% title('accel. from vicon')

% figure
%   %calc rel vel in MP frame
%   vel_rec_rel_GFR = x_rec((1:3)+3+3*pInt,1:N_MP) - x_rec((1:3)+3,1:N_MP);
%   vel_act_rel_GFR = actState(1:N_MP,(1:3)+3+2*pInt) - actState(1:N_MP,(1:3)+3);
%   vel_act_rel_GFR = vel_act_rel_GFR';
% %   if pInt == plotLA
% %       qRotEst = qLankleEst;
% %       qRotAct = qLankle;
% %   %LTIB_CS = quat2rotm(qLankleEst(n,:));
% %   elseif pInt == plotRA
% %       qRotEst = qRankleEst;
% %       qRotAct = qRankle;
% %   %RTIB_CS = quat2rotm(qRankleEst(n,:));
% %   end
% qRotEst = qPelvis;
% qRotAct = qPelvis;
%   %disp(size(vel_rec_rel_GFR))
%   for i = 1:N_MP
%       vel_rec_rel_MPFR(:,i) = quat2rotm(qRotEst(i,:))'*vel_rec_rel_GFR(:,i);
%       vel_act_rel_MPFR(:,i) = quat2rotm(qRotAct(i,:))'*vel_act_rel_GFR(:,i);
%   end
% hold on
% for k=1:3                              % plot results
%   subplot(nDim,1,k)
%   plot((1:N_MP)/fs, vel_rec_rel_MPFR(k,1:N_MP), '-', (1:N_MP)/fs, vel_act_rel_MPFR(k,1:N_MP), '--')
%   ylabel(strcat('state ',num2str(k)))
%   xlabel('time (s)')
%   legend('est vel','tru vel')
% end
% suptitle('rel. vel. in MP frame')
% hold off
%% pos plot
figure
%calc rel pos in MP frame
  pos_est_rel_GFR = estBodyInt(1:N_MP,(1:3)) - estBody.MIDPEL(1:N_MP,1:3);
  pos_act_rel_GFR = actBodyInt(1:N_MP,(1:3)) - actBody.MIDPEL(1:N_MP,1:3);
  %pos_act_rel_GFR = pos_act_rel_GFR';
%   if pInt == plotLA
%       qRotEst = qLankleEst;
%       qRotAct = qLankle;
%   %LTIB_CS = quat2rotm(qLankleEst(n,:));
%   elseif pInt == plotRA
%       qRotEst = qRankleEst;
%       qRotAct = qRankle;
%   %RTIB_CS = quat2rotm(qRankleEst(n,:));
%   end
  %disp(size(vel_rec_rel_GFR))
  for i = 1:N_MP
      pos_est_rel_MPFR(i,:) = quat2rotm(estBody.qRPV(i,:))'*pos_est_rel_GFR(i,:)';
      pos_act_rel_MPFR(i,:) = quat2rotm(actBody.qRPV(i,:))'*pos_act_rel_GFR(i,:)';
  end

for k=1:3                              % plot results
  subplot(nDim,1,k)
  hold on
  plot((1:N_MP)/fs, pos_est_rel_MPFR(1:N_MP,k), '-', (1:N_MP)/fs, pos_act_rel_MPFR(1:N_MP,k), '--')
  ylabel(strcat('state ',num2str(k)))
  xlabel('time (s)')
  legend('est pos','tru pos','act. vel. int.','est. vel. int.')
  set(gca,'fontsize',15)
  box on
end
suptitle(strcat('rel. pos. of',pIntS,'in MP frame'))
hold off
set(gcf, 'Position', [100, 100, plotWidth, plotHeight])
set(gca,'fontsize',15)
box on
filename = strcat(pIntS,'Pos');
formattype = 'png';
saveas(gcf,filename,formattype)

figure
for k=1:3                              % plot results
  subplot(nDim+1,1,k)
  hold on
  plot((1:N_MP)/fs,  pos_act_rel_MPFR(1:N_MP,k) - pos_est_rel_MPFR(1:N_MP,k))
  ylabel(strcat('state ',num2str(k)))
  xlabel('time (s)')
  box on
  set(gca,'fontsize',15)
  %legend('est pos','tru pos','act. vel. int.','est. vel. int.')
end
subplot(nDim+1,1,nDim+1)
ylabel('rms')
pos_res = [pos_act_rel_MPFR(1:N_MP,:) - pos_est_rel_MPFR(1:N_MP,:)]';
pos_res_rms = sqrt(sum(pos_res.^2,1));
plot((1:N_MP)/fs,pos_res_rms)
suptitle(strcat(pIntS,'pos residual (tru - est)'))
hold off
set(gcf, 'Position', [100, 100, plotWidth, plotHeight])
set(gca,'fontsize',15)
box on
filename = strcat(pIntS,'Pos Res');
formattype = 'png';
saveas(gcf,filename,formattype)

disp('rmse = ')
disp(sum(pos_res_rms))
%% additional figure
% plot((1:100)/100, cumtrapz([1:100]/100,[1:100]/100),'.')

% figure
% hold on
% plot3(x_rec(1,1:N_MP),x_rec(2,1:N_MP), x_rec(3,1:N_MP),'o','color','r')
% plot3(x_rec(10,1:N_MP),x_rec(11,1:N_MP), x_rec(12,1:N_MP),'o','color','b')
% plot3(x_rec(19,1:N_MP),x_rec(20,1:N_MP), x_rec(21,1:N_MP),'o','color','c')
% 
% plot3(actState(1:N_MP,1),actState(1:N_MP,2), actState(1:N_MP,3),'.','color','r')
% plot3(actState(1:N_MP,7),actState(1:N_MP,8), actState(1:N_MP,9),'.','color','b')
% plot3(actState(1:N_MP,13),actState(1:N_MP,14), actState(1:N_MP,15),'.', 'color','c')
% view(-90,90)
% title('3d reconstruction')
% hold off
