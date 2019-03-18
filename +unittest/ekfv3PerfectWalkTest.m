%% Test kf_3_kmus_v3 using perfect motion
ideal = struct;
ideal.MIDPEL = [(1:300)' zeros(300,1) 1000*ones(300,1) ] / 1000;
ideal.LFEP = ideal.MIDPEL + ([0 200 0]/1000);
ideal.RFEP = ideal.MIDPEL - ([0 200 0]/1000);
ltheta = [linspace(0, pi/2, 150) linspace(pi/2, 0, 150)]';
rtheta = [linspace(pi/2, 0, 150) linspace(0, pi/2, 150)]';
ideal.LFEO = ideal.LFEP + 0.5*[sin(ltheta) zeros(300,1) -cos(ltheta)];
ideal.RFEO = ideal.RFEP + 0.5*[sin(rtheta) zeros(300,1) -cos(rtheta)];
ideal.LTIO = ideal.LFEO - ([0 0 500]/1000);
ideal.RTIO = ideal.RFEO - ([0 0 500]/1000);
ideal.qRPV = repmat(rotm2quat(eye(3,3)), 300, 1);
ideal.qLTH = zeros(300, 4);
ideal.qRTH = zeros(300, 4);
for i=1:300
    y = [0 1 0];
    z = ideal.RFEP(i,:)-ideal.RFEO(i,:); z = z / norm(z);
    x = cross(y, z);
    ideal.qRTH(i,:) = rotm2quat([x' y' z']);
    
    z2 = ideal.LFEP(i,:)-ideal.LFEO(i,:); z2 = z2 / norm(z2);
    x2 = cross(y, z2);
    ideal.qLTH(i,:) = rotm2quat([x2' y' z2']);
end
ideal.qLSK = repmat(rotm2quat(eye(3,3)), 300, 1);
ideal.qRSK = repmat(rotm2quat(eye(3,3)), 300, 1);

fs = 60;
gfr_vel_MP_act = [gradient(ideal.MIDPEL(:,1), 1/fs), ...
                  gradient(ideal.MIDPEL(:,2), 1/fs), ...
                  gradient(ideal.MIDPEL(:,3), 1/fs)];
gfr_acc_MP_act = [gradient(gfr_vel_MP_act(:,1), 1/fs), ...
                  gradient(gfr_vel_MP_act(:,2), 1/fs), ...
                  gradient(gfr_vel_MP_act(:,3), 1/fs)];
gfr_vel_LA_act = [gradient(ideal.LTIO(:,1), 1/fs), ...
                  gradient(ideal.LTIO(:,2), 1/fs), ...
                  gradient(ideal.LTIO(:,3), 1/fs)];
gfr_acc_LA_act = [gradient(gfr_vel_LA_act(:,1), 1/fs), ...
                  gradient(gfr_vel_LA_act(:,2), 1/fs), ...
                  gradient(gfr_vel_LA_act(:,3), 1/fs)];
gfr_vel_RA_act = [gradient(ideal.RTIO(:,1), 1/fs), ...
                  gradient(ideal.RTIO(:,2), 1/fs), ...
                  gradient(ideal.RTIO(:,3), 1/fs)];
gfr_acc_RA_act = [gradient(gfr_vel_RA_act(:,1), 1/fs), ...
                  gradient(gfr_vel_RA_act(:,2), 1/fs), ...
                  gradient(gfr_vel_RA_act(:,3), 1/fs)];

x0_pos_MP = ideal.MIDPEL(1,:);
x0_pos_LA = ideal.LTIO(1,:);
x0_pos_RA = ideal.RTIO(1,:);
x0_vel_MP = gfr_vel_MP_act(1,:);
x0_vel_LA = gfr_vel_LA_act(1,:);
x0_vel_RA = gfr_vel_RA_act(1,:);
qPelvisEst = ideal.qRPV;
qLankleEst = ideal.qLSK;
qRankleEst = ideal.qRSK;

bIsStatMP = false(300,1);
bIsStatLA = false(300,1);
bIsStatRA = false(300,1);

uwb_mea = struct;
uwb_mea.left_tibia_mid_pelvis = vecnorm((ideal.MIDPEL-ideal.LTIO)', 2, 2);
uwb_mea.mid_pelvis_right_tibia = vecnorm((ideal.MIDPEL-ideal.RTIO)', 2, 2);
uwb_mea.left_tibia_right_tibia = vecnorm((ideal.RTIO-ideal.LTIO)', 2, 2);

d_pelvis = norm(ideal.LFEP(1,:) - ideal.RFEP(1,:));
d_rfemur = norm(ideal.RFEP(1,:) - ideal.RFEO(1,:));
d_lfemur = norm(ideal.LFEP(1,:) - ideal.LFEO(1,:));
d_rtibia = norm(ideal.RFEO(1,:) - ideal.RTIO(1,:));
d_ltibia = norm(ideal.LFEO(1,:) - ideal.LTIO(1,:));

idx = 1:length(ideal.MIDPEL(:,1));
actBody = pelib.grBody('name', 'act', 'posUnit', 'm', 'oriUnit', 'deg', ...
                 'lnSymbol', '-', 'ptSymbol', '.', 'frame', 'world', ...
                 'xyzColor', {'m', 'y', 'c'}, ...
                 'MIDPEL', ideal.MIDPEL(idx,:), ...
                 'LFEP', ideal.LFEP(idx,:), ...
                 'RFEP', ideal.RFEP(idx,:), ...
                 'LFEO', ideal.LFEO(idx,:), ...
                 'RFEO', ideal.RFEO(idx,:), ...
                 'LTIO', ideal.LTIO(idx,:), ...
                 'RTIO', ideal.RTIO(idx,:), ...
                 'qRPV', ideal.qRPV(idx,:), ...
                 'qRTH', ideal.qRTH(idx,:), ...
                 'qLTH', ideal.qLTH(idx,:), ...
                 'qRSK', ideal.qRSK(idx,:), ...
                 'qLSK', ideal.qLSK(idx,:));
actBodyRel = actBody.changeRefFrame('MIDPEL');
actState = [ideal.MIDPEL gfr_vel_MP_act ideal.LTIO gfr_vel_LA_act ...
            ideal.RTIO gfr_vel_RA_act];
resultsIdx = 1; clear results;

for mI = [2]
    for cI = [0 1:8 21:23 51:54 71:78 121:122 124:125 131:132 134:135 ...
              141:144 151:154 201:208 221:223 271:278]
        t0 = cputime;
        
        P = eye(30, 30)*0.5;
        x0 = [x0_pos_MP x0_vel_MP zeros(1,4) ...
              x0_pos_LA x0_vel_LA zeros(1,4) ...
              x0_pos_RA x0_vel_RA zeros(1,4)]';
        v3Options = struct('fs', fs, 'applyCstr', cI, 'applyMeas', mI, ...
            'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5);

        [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.kf_3_kmus_v3(x0, P, ...
            gfr_acc_MP_act, bIsStatMP, qPelvisEst, ...
            gfr_acc_LA_act, bIsStatLA, qLankleEst, ...
            gfr_acc_RA_act, bIsStatRA, qRankleEst, ...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
            v3Options);

        %% --------------------------------------------------------------------
        %  Validation
        estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                               'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'vicon', ...
                               'xyzColor', {'r', 'g', 'b'}, ...
                               'MIDPEL', x_pos_v2(idx, 1:3), ...
                               'LFEP', t_dat_v2.LFEP(idx, :), ...
                               'LFEO', t_dat_v2.LFEO(idx, :), ...
                               'LTIO', x_pos_v2(idx, 11:13), ...
                               'RFEP', t_dat_v2.RFEP(idx, :), ...
                               'RFEO', t_dat_v2.RFEO(idx, :), ...
                               'RTIO', x_pos_v2(idx, 21:23), ...
                               'qRPV', x_pos_v2(idx, 7:10), ...
                               'qLTH', t_dat_v2.qLTH(idx, :), ...
                               'qRTH', t_dat_v2.qRTH(idx, :), ...
                               'qLSK', x_pos_v2(idx, 17:20), ...
                               'qRSK', x_pos_v2(idx, 27:30));

        estBodyRel = estBody.changeRefFrame('MIDPEL');
        results0 = estBodyRel.diffRMSE(actBodyRel);
        results0.name = sprintf('Ideal+M%02d+C%03d', mI, cI);
        results0.runtime = cputime-t0;
        results(resultsIdx) = results0;
        fprintf("Index %3d: %s (%.3f) Running time: %.4f\n", ...
                        resultsIdx, results0.name, results0.posMeanRMSE, cputime-t0);
        
        %% Almost zero position error
        if cI==4 || cI==204
            assert(results0.posMeanRMSE < 0.2)
        else
            assert(results0.posMeanRMSE < 0.02)
        end
        
        %% Check thigh and shank length
        assert(any(abs(vecnorm(estBodyRel.LFEP-estBodyRel.LFEO, 2, 2) - d_lfemur) < 0.01))
        assert(any(abs(vecnorm(estBodyRel.RFEP-estBodyRel.RFEO, 2, 2) - d_rfemur) < 0.01))
        assert(any(abs(vecnorm(estBodyRel.LFEO-estBodyRel.LTIO, 2, 2) - d_ltibia) < 0.01))
        assert(any(abs(vecnorm(estBodyRel.RFEO-estBodyRel.RTIO, 2, 2) - d_rtibia) < 0.01))
        
        %% Check knee angle
        if cI ~= 0
            lkangle = estBodyRel.calcJointAnglesLKnee();
            rkangle = estBodyRel.calcJointAnglesRKnee();
            assert(any(lkangle(:,2) < pi/2))
            assert(any(rkangle(:,2) < pi/2))
        end
        resultsIdx = resultsIdx + 1;
    end
end
results = struct2table(results);

%     updateFigureContents('Animation Freeze 2');
%     grid; view(0, 90); hold on;
%     for i=idx(1):20:idx(end)
%         pelib.viz.plotLowerBody(estBody, i);
%         pelib.viz.plotLowerBody(actBody, i);
%     end
    
% plots
% % Animation
% updateFigureContents('Animation');
% xlabel('x'); ylabel('y'); zlabel('z');
% estBodyLimits = [estBody.xlim() estBody.ylim() estBody.zlim()];
% for i=idx
%     clf; grid;
%     xlim(estBodyLimits(1:2)); 
%     ylim(estBodyLimits(3:4)); 
%     zlim(estBodyLimits(5:6));  
%     view(-126, 26);
%     pelib.viz.plotLowerBody(estBody, i);
%     pause(1/1000);
% end
% 
% updateFigureContents('Animation');
% xlabel('x'); ylabel('y'); zlabel('z');
% actBodyLimits = [actBody.xlim() actBody.ylim() actBody.zlim()];
% for i=idx
%     clf; grid; 
%     xlim(actBodyLimits(1:2)); 
%     ylim(actBodyLimits(3:4)); 
%     zlim(actBodyLimits(5:6));  
%     view(-126, 26);
%     pelib.viz.plotLowerBody(actBody, i);
%     pause(1/1000);
% end
% 
% % Static Plots
% updateFigureContents('Position');
% pelib.viz.plotPosition({estBody, actBody}, {'MIDPEL', 'LTIO', 'RTIO'});