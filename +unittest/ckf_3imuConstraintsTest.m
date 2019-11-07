%% ckf_3imus hinge joint constraints test
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
      
%% Ideal case straight
d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
fs = 60;
mIList = [2];
cIList = [355];
      
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
cIList = [355];
      
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
            if cI==5 || cI==75 || cI == 205 || cI == 275
                threshold = 0.2;
            else
                threshold = 0.01;
            end
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < threshold))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < threshold))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < threshold))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < threshold))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < threshold))
        end
    end
end

%% Knee inequality check
d_pelvis = 0.5;
d_lfemur = 0.6; d_rfemur = 0.6;
d_ltibia = 0.5; d_rtibia = 0.5;
fs = 60;
mIList = [0];
cIList = [355];
      
xNextList = [struct('MIDPEL', [0.3 0 1.0], ...
                    'LFEP', [0.3 0.25 1.0], 'RFEP', [0.3 -0.25 1.0], ...
                    'LFEO', [0.0 0.25 0.5], 'RFEO', [0.0 -0.25 0.5], ...
                    'LTIO', [0 0.25 0], 'RTIO', [0 -0.25 0], ...
                    'LKAngle', 0, 'RKAngle', 0, ...
                    'LKAngleAct', 0, 'RKAngleAct', 0), ...
             struct('MIDPEL', [0 0 1.1], ...
                    'LFEP', [0 0.25 1.1], 'RFEP', [0 -0.25 1.1], ...
                    'LFEO', [0 0.25 0.5], 'RFEO', [0 -0.25 0.5], ...
                    'LTIO', [0 0.25 1.0], 'RTIO', [0 -0.25 1.0], ...
                    'LKAngle', 160, 'RKAngle', 160, ...
                    'LKAngleAct', 160, 'RKAngleAct', 160), ...
             ];
for xNext = xNextList
    for mI = mIList
        for cI = cIList
            threshold = 0.01;
            estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                                        d_ltibia, d_rtibia, fs, mI, cI);
            assert(any(abs(vecnorm(estBody.LFEP-estBody.RFEP, 2, 2) - d_pelvis) < threshold))
            assert(any(abs(vecnorm(estBody.LFEP-estBody.LFEO, 2, 2) - d_lfemur) < threshold))
            assert(any(abs(vecnorm(estBody.RFEP-estBody.RFEO, 2, 2) - d_rfemur) < threshold))
            assert(any(abs(vecnorm(estBody.LFEO-estBody.LTIO, 2, 2) - d_ltibia) < threshold))
            assert(any(abs(vecnorm(estBody.RFEO-estBody.RTIO, 2, 2) - d_rtibia) < threshold))
            
%             lkangle0 = calcKneeAngle(xNext.LFEP, xNext.LFEO, xNext.LTIO);
%             rkangle0 = calcKneeAngle(xNext.RFEP, xNext.RFEO, xNext.RTIO);
            lkangle = estBody.calcJointAnglesLKnee()*180/pi;
            rkangle = estBody.calcJointAnglesRKnee()*180/pi;
            assert(abs(lkangle(:,2)-xNext.LKAngleAct) < 1e-2)
            assert(abs(rkangle(:,2)-xNext.RKAngleAct) < 1e-2)
        end
    end
end

function theta = calcKneeAngle(Hip, Knee, Ankle)
    x = Hip - Knee;
    y = Ankle - Knee;
    theta = acosd(dot(x, y) / (norm(x)*norm(y)));
end

function estBody = runConstraintTest(xNext, d_pelvis, d_lfemur, d_rfemur, ...
                    d_ltibia, d_rtibia, fs, mI, cI)
    v3Options = struct('fs', fs, 'applyCstr', cI, 'applyMeas', mI, ...
            'sigma2QAccMP', 0.5^2, 'sigma2QAccLA', 0.5^2, ...
            'sigma2QAccRA', 0.5^2 );
        
    x0 = [0 0 0 xNext.MIDPEL*fs ...
          0 0 0 xNext.LTIO*fs ...
          0 0 0 xNext.RTIO*fs ]';
    P = eye(18, 18)*0.5;
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

    [ x_pri_v2, x_pos_v2, t_dat_v2 ] = pelib.est.ckf_3imus(x0, P, ...
            gfr_acc_MP_act, bIsStatMP, qPelvisEst, ...
            gfr_acc_LA_act, bIsStatLA, qLankleEst, ...
            gfr_acc_RA_act, bIsStatRA, qRankleEst, ...
            d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
            v3Options);

    idx = 1:1;
    estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
                   'xyzColor', {'r', 'g', 'b'}, 'fs', fs, ...
                   'MIDPEL', x_pos_v2(idx, 1:3), ...
                   'LFEP', t_dat_v2.LFEP(idx, :), ...
                   'LFEO', t_dat_v2.LFEO(idx, :), ...
                   'LTIO', x_pos_v2(idx, 7:9), ...
                   'RFEP', t_dat_v2.RFEP(idx, :), ...
                   'RFEO', t_dat_v2.RFEO(idx, :), ...
                   'RTIO', x_pos_v2(idx, 13:15), ...
                   'qRPV', qPelvisEst(idx, :), ...
                   'qLTH', t_dat_v2.qLTH(idx, :), ...
                   'qRTH', t_dat_v2.qRTH(idx, :), ...
                   'qLSK', qLankleEst(idx, :), ...
                   'qRSK', qRankleEst(idx, :));
end