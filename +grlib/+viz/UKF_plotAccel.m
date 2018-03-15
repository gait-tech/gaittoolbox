function UKF_plotAccel(pIntS,tru_Acc_GFR, est_Acc_GFR, N_MP, fs)
plotHeight = 800;
plotWidth = 1000;
figure
for k=1:3                              % plot results
    subplot(3,1,k)
    hold on
    plot((1:N_MP)/fs, est_Acc_GFR(1:N_MP,k), '-', (1:N_MP)/fs, tru_Acc_GFR(1:N_MP,k), '--')
    ylabel(strcat('state ',num2str(k)))
    xlabel('time (s)')
    legend('est acc.','tru acc.')
    box on
    set(gca,'fontsize',15)
end
suptitle(strcat('accel of ',pIntS,' in GFR'))
hold off
set(gcf, 'Position', [100, 100, plotWidth, plotHeight])
set(gca,'fontsize',15)
filename = strcat(pIntS,'Acc');
formattype = 'png';
saveas(gcf,filename,formattype)

figure
box on
for k=1:3                              % plot results
    subplot(3+1,1,k)
    hold on
    plot((1:N_MP)/fs,  tru_Acc_GFR(1:N_MP,k) - est_Acc_GFR(1:N_MP,k))
    ylabel(strcat('state ',num2str(k)))
    xlabel('time (s)')
    box on
    set(gca,'fontsize',15)
end
subplot(3+1,1,3+1)
ylabel('rms')
pos_res = [tru_Acc_GFR(1:N_MP,:) - est_Acc_GFR(1:N_MP,:)]';
pos_res_rms = sqrt(sum(pos_res.^2,1));
plot((1:N_MP)/fs,pos_res_rms)
suptitle(strcat(pIntS,'acc. residual (tru - est)'))
hold off
set(gcf, 'Position', [100, 100, plotWidth, plotHeight])
set(gca,'fontsize',15)
box on
filename = strcat(pIntS,'Acc Res');
formattype = 'png';
saveas(gcf,filename,formattype)