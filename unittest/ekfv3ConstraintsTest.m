%% kf_3_kmus hinge joint constraints test
% Author: Luke Sy
%
% play with the hinge joint constraint and understand its behavior
% 
% Findings:
% 1. The update ('d') is only based on position error. If P=I, only the
% position states will be updated depending on P (higher value means less 
% trust worthy. The velocity states will only be updated due to covariance 
% between position and velocity.
% 2. The update does not change the knee angle (alpha). Hence, if position
% error causes knee angle error, the error is carried over and not fixed by
% this constraint update.
path(path, genpath('..'));
      
%% Ideal case straight
d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
fs = 60;
mIList = [2];
cIList = [1:8 21:23 51:54 71:78 121:122 124:125 131:132 134:135 ...
          141:144 151:154 201:208 221:223 271:278];
      
xNext = struct('MIDPEL', [0 0 1.1], ...
               'LFEO', [0 0.25 0.5], 'RFEO', [0 -0.25 0.5], ...
               'LTIO', [0 0.25 0], 'RTIO', [0 -0.25 0] );
for mI = mIList
    for cI = cIList
        estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                    d_ltibia, d_rtibia, fs, mI, cI);
        assert(any(estBody.MIDPEL == xNext.MIDPEL))
        assert(any(estBody.LTIO == xNext.LTIO))
        assert(any(estBody.RTIO == xNext.RTIO))
    end
end

%% Longer or shorter limbs
d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
fs = 60;
mIList = [2];
cIList = [1:8 21:23 51:54 71:78 121:122 124:125 131:132 134:135 ...
          141:144 151:154 201:208 221:223 271:278];
      
xNextList = [struct('MIDPEL', [0 0 1.2], ...
                    'LFEO', [0 0.3 0.7], 'RFEO', [0 -0.3 0.5], ...
                    'LTIO', [0 0.3 0.2], 'RTIO', [0 -0.3 0] ), ...
             struct('MIDPEL', [0 0 1.1], ...
                    'LFEO', [d_lfemur*sin(pi/4) 0.25 1.1-d_lfemur*cos(pi/4)], ...
                    'RFEO', [d_rfemur*sin(pi/4) 0.25 1.1-d_rfemur*cos(pi/4)], ...
                    'LTIO', [d_lfemur*sin(pi/4) 0.25 1.1-d_lfemur*cos(pi/4)-d_ltibia], ...
                    'RTIO', [d_lfemur*sin(pi/4) -0.25 1.1-d_rfemur*cos(pi/4)-d_rtibia] ) ...
             ];
for xNext = xNextList
    for mI = mIList
        for cI = cIList
            estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                        d_ltibia, d_rtibia, fs, mI, cI);
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < 0.01))
        end
    end
end

%% Knee inequality check
d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
fs = 60;
mIList = [2];
cIList = [1:8 21:23 51:54 71:78 121:122 124:125 131:132 134:135 ...
          141:144 151:154 201:208 221:223 271:278];
      
xNextList = [struct('MIDPEL', [0.3 0 1.0], ...
                    'LFEO', [0 0.25 0.5], 'RFEO', [0 -0.25 0.5], ...
                    'LTIO', [0 0.25 0], 'RTIO', [0 -0.25 0], ...
                    'LKAngle', 0, 'RKAngle', 0, ...
                    'LKAngleAct', -30.9638, 'RKAngleAct', -30.9638), ...
             struct('MIDPEL', [0 0 1.1], ...
                    'LFEO', [0 0.25 0.5], 'RFEO', [0 -0.25 0.5], ...
                    'LTIO', [0 0.25 1.0], 'RTIO', [0 -0.25 1.0], ...
                    'LKAngle', 160, 'RKAngle', 160, ...
                    'LKAngleAct', 180, 'RKAngleAct', 180), ...
             ];
for xNext = xNextList
    for mI = mIList
        for cI = [1:8 21:23 51:54 71:78 121:122 124:125 141:144]
            estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                        d_ltibia, d_rtibia, fs, mI, cI);
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < 0.01))
            
            lkangle = estBody.calcJointAnglesLKnee()*180/pi;
            rkangle = estBody.calcJointAnglesRKnee()*180/pi;
            assert(abs(lkangle(:,2)-xNext.LKAngleAct) < 1e-2)
            assert(abs(rkangle(:,2)-xNext.RKAngleAct) < 1e-2)
        end
        
        for cI = [131 134]
            estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                        d_ltibia, d_rtibia, fs, mI, cI);
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < 0.01))
            
            lkangle = estBody.calcJointAnglesLKnee()*180/pi;
            rkangle = estBody.calcJointAnglesRKnee()*180/pi;
            assert(any(lkangle(:,2) < 180 && lkangle(:,2) > 0))
            assert(any(rkangle(:,2) < 180 && rkangle(:,2) > 0))
        end
        
        for cI = [132 135 151:154 201:208 221:223 271:278]
            estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                        d_ltibia, d_rtibia, fs, mI, cI);
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < 0.01))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < 0.01))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < 0.01))
            
            lkangle = estBody.calcJointAnglesLKnee()*180/pi;
            rkangle = estBody.calcJointAnglesRKnee()*180/pi;
            assert(abs(lkangle(:,2)-xNext.LKAngle) < 0.5)
            assert(abs(rkangle(:,2)-xNext.RKAngle) < 0.5)
        end
    end
end
% 
% % Ideal case bent
% % xhat = [0 0 0 0.1 0.1 0 ...
% %         0.6 0.25 -0.5 0.1 -0.1 0 ...
% %         0.6 -0.25 -0.5 0.1 -0.1 0]';
% % Non-ideal case bent, longer/shorter length
% % xhat = [0 0 0 0.1 0.1 0 ...
% %         0.7 0.25 -0.7 0.1 -0.1 0 ...
% %         0.7 -0.25 -0.4 0.1 -0.1 0]';
% 
% % Ideal case slightly bent
% % xhat = [0 0 0 0.1 0.1 0 ...
% %         0.6*sin(pi/4) 0.25 0.6*cos(pi/4)-0.5 0.1 -0.1 0 ...
% %         0.6*sin(pi/4) -0.25 0.6*cos(pi/4)-0.5 0.1 -0.1 0]';
% % Ideal case slightly bent, knee not hinge joint
% % xhat = [0 0 0 0.1 0.1 0 ...
% %         0.6*sin(pi/4) 0.5 0.6*cos(pi/4)-0.5 0.1 -0.1 0 ...
% %         0.6*sin(pi/4) -0.5 0.6*cos(pi/4)-0.5 0.1 -0.1 0]';
% 
% LKNE = xhat(7:9,1) + d_ltibia*LTIB_CS(:,3);
% RKNE = xhat(13:15,1) + d_rtibia*RTIB_CS(:,3);
% LFEM_z = xhat(1:3,1)+d_pelvis/2*PELV_CS(:,2)-LKNE;
% RFEM_z = xhat(1:3,1)-d_pelvis/2*PELV_CS(:,2)-RKNE;
% alpha_lk = acos(dot(LFEM_z, LTIB_CS(:,3))/(norm(LFEM_z)*norm(LTIB_CS(:,3))));
% alpha_rk = acos(dot(RFEM_z, RTIB_CS(:,3))/(norm(RFEM_z)*norm(RTIB_CS(:,3))));
% 
% d_k = [ (d_pelvis/2*PELV_CS(:,2) ...
%          -d_lfemur*cos(alpha_lk)*LTIB_CS(:,3) ...
%          +d_lfemur*sin(alpha_lk)*LTIB_CS(:,1) ...
%          -d_ltibia*LTIB_CS(:,3)) ; ...
%         (-d_pelvis/2*PELV_CS(:,2)+ ...
%          -d_rfemur*cos(alpha_rk)*RTIB_CS(:,3) ...
%          +d_rfemur*sin(alpha_rk)*RTIB_CS(:,1) ...
%          -d_rtibia*RTIB_CS(:,3)) ];
% 
% res = d_k-D*xhat;
% Kk = P*D'*(D*P*D')^(-1);
% dx = Kk*(res);
% xhat = xhat + dx;
% 
% LKNE = xhat(7:9,1) + d_ltibia*LTIB_CS(:,3);
% RKNE = xhat(13:15,1) + d_rtibia*RTIB_CS(:,3);
% LHIP = xhat(1:3,1)+d_pelvis/2*PELV_CS(:,2);
% RHIP = xhat(1:3,1)-d_pelvis/2*PELV_CS(:,2);
% LFEM_z = LHIP-LKNE;
% RFEM_z = RHIP-RKNE;
% alpha_lk = acos(dot(LFEM_z, LTIB_CS(:,3))/(norm(LFEM_z)*norm(LTIB_CS(:,3))));
% alpha_rk = acos(dot(RFEM_z, RTIB_CS(:,3))/(norm(RFEM_z)*norm(RTIB_CS(:,3))));

function estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                    d_ltibia, d_rtibia, fs, mI, cI)
    v3Options = struct('fs', fs, 'applyCstr', cI, 'applyMeas', mI, ...
            'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5);
        
    x0 = [0 0 0 xNext.MIDPEL*fs 0 0 0 0 ...
          0 0 0 xNext.LTIO*fs 0 0 0 0 ...
          0 0 0 xNext.RTIO*fs 0 0 0 0 ]';
    P = eye(30, 30)*0.5;
    gfr_acc_MP_act = zeros(1, 3);
    gfr_acc_LA_act = zeros(1, 3);
    gfr_acc_RA_act = zeros(1, 3);
    qPelvisEst = rotm2quat(eye(3,3));
    y = [0 1 0];
    z = xNext.LFEO - xNext.LTIO; z = z/norm(z);
    x = cross(y, z);
    qLankleEst = rotm2quat([x' y' z']);
    z = xNext.RFEO - xNext.RTIO; z = z/norm(z);
    x = cross(y, z);
    qRankleEst = rotm2quat([x' y' z']);
    bIsStatMP = zeros(1, 1); bIsStatLA = zeros(1, 1); bIsStatRA = zeros(1, 1);
    uwb_mea = struct;
    uwb_mea.left_tibia_mid_pelvis = zeros(1, 1);
    uwb_mea.mid_pelvis_right_tibia = zeros(1, 1);
    uwb_mea.left_tibia_right_tibia = zeros(1, 1);
    
    [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.kf_3_kmus_v3(x0, P, ...
            gfr_acc_MP_act, bIsStatMP, qPelvisEst, ...
            gfr_acc_LA_act, bIsStatLA, qLankleEst, ...
            gfr_acc_RA_act, bIsStatRA, qRankleEst, ...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
            v3Options);

    idx = 1:1;
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

end