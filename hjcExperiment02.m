% check effect of putting acc bias in state

gfr_acc_MP2 = gfr_acc_MP(1:end-1,:);
gfr_acc_MP_nb = gfr_acc_MP2 - ...
    quatrotate(quatconj(qPelvisAct(2:end,:)), x_pos_v2(:,19:21));

gfr_acc_LA2 = gfr_acc_LA(1:end-1,:);
gfr_acc_LA_nb = gfr_acc_LA2 - ...
    quatrotate(quatconj(qLankleAct(2:end,:)), x_pos_v2(:,22:24));

gfr_acc_RA2 = gfr_acc_RA(1:end-1,:);
gfr_acc_RA_nb = gfr_acc_RA2 - ...
    quatrotate(quatconj(qRankleAct(2:end,:)), x_pos_v2(:,25:27));

[ pelib.rmse(gfr_acc_MP_act-gfr_acc_MP_nb) ...
  pelib.rmse(gfr_acc_MP_act-gfr_acc_MP2); ...
  pelib.rmse(gfr_acc_LA_act-gfr_acc_LA_nb) ...
  pelib.rmse(gfr_acc_LA_act-gfr_acc_LA2); ...
  pelib.rmse(gfr_acc_RA_act-gfr_acc_RA_nb) ...
  pelib.rmse(gfr_acc_RA_act-gfr_acc_RA2) ]

pelib.viz.plotXYZ(gfr_acc_MP_act, gfr_acc_MP_nb, gfr_acc_MP)
pelib.viz.plotXYZ(gfr_acc_MP_act, gfr_acc_MP)
pelib.viz.plotXYZ(gfr_acc_MP_act, gfr_acc_MP_nb)
pelib.viz.plotXYZ(gfr_acc_MP, gfr_acc_MP_nb)


