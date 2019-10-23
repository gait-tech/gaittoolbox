function [ xhat_pri, xhat_con, debug_dat ] = ckf_3imus(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, options)
% Constrained Kalman Filter for performing sensor fusion on the trajectory of
% three KMUs presumably worn on the body in the following configuration: mid
% pelvis, left ankle, right ankle
% In this state space model, the position and velocity of each kinematic
% measurement unit (KMU) is estimated in 3D space by combining the
% information from each KMU in a kalman filter. NOTE: pay special attention 
% to units:;
% position (meters)
% velocity (m/s)
% acceleration (m/2^2)
%
%   :param x0: the initial state in the GFR
%   :param gfrAccMP: the acceleration of the mid-pelvis in the GFR
%   :param gfrAccLA: the acceleration of the left ankle in the GFR
%   :param gfrAccRA: the acceleration of the right ankle in the GFR
%   :param bIsStatMP: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param bIsStatLA: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param bIsStatRA: a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   :param qMP:		mid  pelvis orientation in the GFR (quaternion)
%   :param qLA:		left  ankle orientation in the GFR (quaternion)
%   :param qRA:		right ankle orientation in the GFR (quaternion)
%   :param dPelvis:	pelvis width
%   :param dRFemur:	right femur length
%   :param dLFemur:	left femur length
%   :param dRTibia:	right tibia length
%   :param dLTibia:	left tibia length
%   :param options:	struct containing the following settings:
%       applyMeas - turn on/off zero velocity update. boolean
%           X=0: normal measurement update
%           X=1: add 3 point dist
%           X=2: add calc knee angle from uwb dist
%           X=3: add calc knee angle from 3 point dist
%           X01: standard zupt
%           X02: standard zupt + floor assumption
%           X03: standard zupt + floor assumption + reset at both foot
%           X04: standard zupt + floor assumption + series reset at both foot
%           X05: standard zupt + floor assumption + zero first floor foot
%           X06: standard zupt + floor assumption + pelvis z pos assumption
%           X07: standard zupt + floor assumption + pelvis z pos assumption + 
%                zero first floor foot
%           X08: standard zupt + floor assumption + reset first floor foot
%           X21: standard zupt + floor assumption + estimate projection (W=P^-1) assuming perfect orientation
%           X31: standard zupt + floor assumption + recenter at left foot step
%                + one big cov update
%           X32: standard zupt + floor assumption + limit pos max covariance
%                + one big cov update
%           X33: standard zupt + floor assumption + recenter at left foot step
%                + limit pos max covariance + one big cov update
%           X41: standard zupt + floor assumption + recenter at left foot step
%                + series cov update
%           X42: standard zupt + floor assumption + limit pos max covariance
%                + series cov update
%           X43: standard zupt + floor assumption + recenter at left foot step
%                + limit pos max covariance + series cov update
%           X50-X57: standard zupt + floor assumption during step & below floor + 
%           X60-X67: standard zupt + floor assumption during step + 
%           X70-X77: standard zupt + floor assumption during step 
%                    + limit pos max covariance with one big cov update
%           X80-X87: standard zupt + floor assumption during step + below floor zupt
%                    + limit pos max covariance with one big cov update
%                    1st bit (+1): pelvis = ankle speed
%                    2nd bit (+2): pelvis = ankle XY pos
%                    3rd bit (+4): pelvis = initial pelvis z pos
%       applyCstr - turn on/off constraints.
%           001: estimate projection (W=P^-1) assuming perfect orientation
%           002: estimate projection (W=I) assuming perfect orientation
%           003: least squares estimate w/ full confidence on pelvis (P=0 at pelvis) + no P update
%           004: maximum probability estimate w/ force equal foot covariance + no P update
%           005: maximum probability estimate if P is well conditioned + P update
%           006: soft maximum probability estimate + P update
%           007: maximum probability estimate of constraint subset + P update
%           008: maximum probability estimate + soft P update
%           011: fmincon (interior point) linear hjc in world frame, W = P^-1
%           012: fmincon (sqp) linear hjc in world frame, W = P^-1
%           013: fmincon (active set) linear hjc in world frame, W = P^-1
%           014: fmincon (interior point) linear hjc in world frame, W = I
%           015: fmincon (sqp) linear hjc in world frame, W = I
%           016: fmincon (active set) linear hjc in world frame, W = I
%           021: PELVIS frame constraint update, maximum probability estimate + no P update
%           022: PELVIS frame constraint update, least squares estimate + no P update
%           023: maximum probability estimate of constraint subset + P update
%           031: fmincon (interior point) linear hjc in pelvis frame, W = P^-1
%           032: fmincon (sqp) linear hjc in pelvis frame, W = P^-1
%           033: fmincon (active set) linear hjc in pelvis frame, W = P^-1
%           034: fmincon (interior point) linear hjc in pelvis frame, W = I
%           035: fmincon (sqp) linear hjc in pelvis frame, W = I
%           036: fmincon (active set) linear hjc in pelvis frame, W = I
%           051: estimate projection (W=P^-1) assuming perfect orientation 
%                + MP/LA/RA zpos = floor zpos
%           052: estimate projection (W=I) assuming perfect orientation
%                + MP/LA/RA zpos = floor zpos
%           053: soft maximum probability estimate 
%			     + MP/LA/RA zpos = floor zpos + P update
%           054: maximum probability estimate of constraint subset 
%                + MP/LA/RA zpos = floor zpos + P update
%           071: estimate projection (W=P^-1) assuming perfect orientation
%                + lowest point = floor
%           072: estimate projection (W=I) assuming perfect orientation
%                + lowest point = floor
%           073: least squares estimate w/ full confidence on pelvis (P=0 at pelvis) + no P update
%                + lowest point = floor
%           074: maximum probability estimate w/ force equal foot covariance + no P update
%                + lowest point = floor
%           075: maximum probability estimate if P is well conditioned + P update
%                + lowest point = floor
%           076: soft maximum probability estimate + P update
%                + lowest point = floor
%           077: maximum probability estimate of constraint subset + P update
%                + lowest point = floor
%           078: maximum probability estimate + soft P update
%                + lowest point = floor
%           101: nonlinear fmincon interior-point (W=P^-1) + knee lock + no P update
%           102: nonlinear fmincon sqp (W=P^-1) + knee lock + no P update
%           103: nonlinear fmincon active set (W=P^-1) + knee lock + no P update
%           104: nonlinear fmincon interior-point (W=I) + knee lock  + no P update
%           105: nonlinear fmincon sqp (W=I) + knee lock + no P update
%           106: nonlinear fmincon active set (W=I) + knee lock + no P update
%           111: nonlinear fmincon interior-point (W=P^-1) + knee lock and ineq + no P update
%           112: nonlinear fmincon sqp (W=P^-1) + knee lock and ineq + no P update
%           113: nonlinear fmincon active set (W=P^-1) + knee lock and ineq + no P update
%           114: nonlinear fmincon interior-point (W=I) + knee lock and ineq + no P update
%           115: nonlinear fmincon sqp (W=I) + knee lock and ineq + no P update
%           116: nonlinear fmincon active set (W=I) + knee lock and ineq + no P update
%           121: nonlinear fmincon interior-point (W=P^-1) + x=pos only 
%                + adj knee + no P update
%           122: nonlinear fmincon sqp (W=P^-1) + x=pos only
%                + adj knee + no P update
%           123: nonlinear fmincon active set (W=P^-1) + x=pos only 
%                + adj knee + no P update
%           124: nonlinear fmincon interior-point (W=I) + x=pos only
%                + adj knee + no P update
%           125: nonlinear fmincon sqp (W=I) + x=pos only 
%                + adj knee + no P update
%           126: nonlinear fmincon active set (W=I) + x=pos only 
%                + adj knee + no P update
%           131: nonlinear fmincon interior-point (W=P^-1) + x=pos only 
%                + adj knee + knee ineq + no P update
%           132: nonlinear fmincon sqp (W=P^-1) + x=pos only 
%                + adj knee + knee ineq + no P update
%           133: nonlinear fmincon active set (W=P^-1) + x=pos only 
%                + adj knee + knee ineq + no P update
%           134: nonlinear fmincon interior-point (W=I) + x=pos only 
%                + adj knee + knee ineq + no P update
%           135: nonlinear fmincon sqp (W=I) + x=pos only 
%                + adj knee + knee ineq + no P update
%           136: nonlinear fmincon active set (W=I) + x=pos only 
%                + adj knee + knee ineq + no P update
%           141: smoothly constraint kf (W=P^-1, early stop) + adj knee + no P update
%           142: smoothly constraint kf (W=I, early stop) + adj knee + no P update
%           143: smoothly constraint kf (W=P^-1, use maxIter) + adj knee + no P update
%           144: smoothly constraint kf (W=I, use maxIter) + adj knee + no P update
%           151: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update
%           152: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update
%           153: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update
%           154: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update
%           155: smoothly constraint kf (W=pelvis frame, early stop)
%                + adj knee + knee ineq + no P update
%           156: smoothly constraint kf (W=pelvis frame, use maxIter)
%                + adj knee + knee ineq + no P update
%           157: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%           158: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P update
%           161: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq
%           162: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq
%           163: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq
%           164: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq
%           165: smoothly constraint kf (W=pelvis frame, early stop)
%                + adj knee + knee ineq
%           166: smoothly constraint kf (W=pelvis frame, use maxIter)
%                + adj knee + knee ineq
%           167: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq
%           168: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq
%           171: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           172: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           173: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           174: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           175: smoothly constraint kf (W=pelvis frame, early stop)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           176: smoothly constraint kf (W=pelvis frame, use maxIter)
%                + adj knee + knee ineq + no P update + lowest point = floor
%           177: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%                + lowest point = floor
%           178: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P update
%                + lowest point = floor
%           201: estimate projection (W=P^-1) assuming perfect orientation
%                + knee angle inequality constraint
%           202: estimate projection (W=I) assuming perfect orientation
%                + knee angle inequality constraint
%           203: least squares estimate w/ full confidence on pelvis (P=0 at pelvis) 
%                + no P update + knee angle inequality constraint
%           204: maximum probability estimate w/ force equal foot covariance 
%                + no P update + knee angle inequality constraint
%           205: maximum probability estimate if P is well conditioned 
%                + P update + knee angle inequality constraint
%           206: soft maximum probability estimate + P update
%                + knee angle inequality constraint
%           207: maximum probability estimate of constraint subset 
%                + P update + knee angle inequality constraint
%           208: maximum probability estimate + soft P update
%           221: PELVIS frame constraint update, maximum probability estimate 
%                + no P update + knee angle inequality constraint
%           222: PELVIS frame constraint update, least squares estimate 
%                + no P update + knee angle inequality constraint
%           223: maximum probability estimate of constraint subset 
%                + P update + knee angle inequality constraint
%           271: estimate projection (W=P^-1) assuming perfect orientation
%                + lowest point = floor + knee angle inequality constraint
%           272: estimate projection (W=I) assuming perfect orientation
%                + lowest point = floor + knee angle inequality constraint
%           273: least squares estimate w/ full confidence on pelvis (P=0 at pelvis) + no P update
%                + lowest point = floor + knee angle inequality constraint
%           274: maximum probability estimate w/ force equal foot covariance + no P update
%                + lowest point = floor + knee angle inequality constraint
%           275: maximum probability estimate if P is well conditioned + P update
%                + lowest point = floor + knee angle inequality constraint
%           276: soft maximum probability estimate + P update
%                + lowest point = floor + knee angle inequality constraint
%           277: maximum probability estimate of constraint subset + P update
%                + lowest point = floor + knee angle inequality constraint
%           278: maximum probability estimate + soft P update
%                + lowest point = floor + knee angle inequality constraint
%           WX1: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update
%           WX2: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update
%           WX3: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update
%           WX4: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update
%           WX5: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%           WX6: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P updates
%                W=3: Use all 30 states
%                W=4: Use only 12 states (pelv, lank, rank position)
%                W=5: uwb meas + use all 30 states
%                X=1: foot step (z pos) = floor (step detect only)
%                X=2: foot step (z pos) = floor (step detect only)
%                     + static foot step (x,y pos)
%                X=3: foot step (z pos) = floor (step detect only)
%                     + zero sigma for foot x,y pos
%                X=5: no foot related locking
%                X=7: foot step (z pos) = floor (step detect and < floorZ)
%                X=8: foot step (z pos) = floor (step detect and < floorZ)
%                     + static foot step (x,y pos)
%                X=9: foot step (z pos) = floor (step detect and < floorZ)
%                     + zero sigma for foot x,y pos
% 
% .. Author: - Luke Wicent Sy, Michael Del Rosario
    fOpt = struct('fs', 60, 'applyMeas', 76, 'applyCstr', 355, ...
        'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 1e3, 'sigmaQOriLA', 1e3, 'sigmaQOriRA', 1e3, ...
        'sigmaROriMP', 1e-3, 'sigmaROriLA', 1e-3, 'sigmaROriRA', 1e-3, ...
        'sigmaRPosLA', 1e-2, 'sigmaRPosRA', 1e-2, 'sigmaRPosZMP', 1e-1, ...
        'sigmaRPosMPLimit', 1e2, 'sigmaRPosLALimit', 1e2, ...
        'sigmaRPosRALimit', 1e2, 'sigmaRPosLARecenter', 1e-3, ...
        'sigmaRVelMPLARA', 1e1, 'sigmaRPosMPLARA', 1e1, ...
        'sigmaCPos', 1e-2, ...
        'sigmaUwbMPLA', 1e1, 'sigmaUwbMPRA', 1e1, 'sigmaUwbLARA', 1e1, ...
        'sigmaUwbLLeg', 1e0, 'sigmaUwbRLeg', 1e0, ...
        'sigmaZuptMP', 1e-1, 'sigmaZuptLA', 1e-1, 'sigmaZuptRA', 1e-1, ...
        'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
        'alphaRKmin', 0, 'alphaRKmax', pi*8/9, ...
        'optimOptimalityTolerance', 1e-2, ...
        'optimConstraintTolerance', 1e-2, ...
        'optimMaxFunctionEvaluations', 1500, 'optimUseParallel', false, ...
        'sckfAlpha', 0.1, 'sckfThreshold', 100, 'sckfMaxIter', 500);
    
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    idxVelLA = 14:16; % column idx corresponding to the left ankle velocity
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    idxVelRA = 24:26; % column idx corresponding to the right ankle velocity
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation
    idxMOriMP = 1:4;
    idxMOriLA = 5:8;
    idxMOriRA = 9:12;
    
    nStates = 30;
  
    if fOpt.applyAccBias
        idxAccBiasMP = 31:33; % column idx corresponding to the mid-pelvis acc bias
        idxAccBiasLA = 34:36; % column idx corresponding to the left ankle acc bias
        idxAccBiasRA = 37:39; % column idx corresponding to the right ankle acc bias
        nStates = 39;
    end
    
    % initialise state vector (must be column)
    validateattributes(x0, {'numeric'}, ...
                       {'2d', 'ncols', 1, 'nrows', nStates});
    % demo start on floor start
%     lowestpoint = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
%     x0(idxPosMP(3)) = x0(idxPosMP(3)) - lowestpoint;
%     x0(idxPosLA(3)) = x0(idxPosLA(3)) - lowestpoint;
%     x0(idxPosRA(3)) = x0(idxPosRA(3)) - lowestpoint;
    % demo start on floor end
    x_tilde = x0;

    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt*dt;        % local variable for readability
    I_N = eye(nStates);
    
    % state transition matrix encodes the relationship between previous state
    % estimate and current state estimate
    F = eye(nStates,nStates);
    % x = x(t-1) + v(t-1)*dt + 0.5*a(t)*dt^2
    % x = A*x + B*u
    F(idxPosMP,idxVelMP) = dt.*eye(3); % mid pelvis
    F(idxPosLA,idxVelLA) = dt.*eye(3); % left ankle
    F(idxPosRA,idxVelRA) = dt.*eye(3); % right ankle

    if fOpt.applyAccBias
        F(idxPosMP, idxAccBiasMP) = -dt2.*eye(3);
        F(idxPosLA, idxAccBiasLA) = -dt2.*eye(3);
        F(idxPosRA, idxAccBiasRA) = -dt2.*eye(3);
        F(idxVelMP, idxAccBiasMP) = -dt.*eye(3);
        F(idxVelLA, idxAccBiasLA) = -dt.*eye(3);
        F(idxVelRA, idxAccBiasRA) = -dt.*eye(3);
    end
    
    G = zeros(nStates,9);
    G(idxPosMP, 1:3) = dt2.*eye(3);
    G(idxVelMP, 1:3) = dt .*eye(3);
    G(idxPosLA, 4:6) = dt2.*eye(3);
    G(idxVelLA, 4:6) = dt .*eye(3);
    G(idxPosRA, 7:9) = dt2.*eye(3);
    G(idxVelRA, 7:9) = dt .*eye(3);

    % Initialise process noise covariance
    Q = diag(repelem([(fOpt.sigmaQAccMP)^2 (fOpt.sigmaQAccLA)^2 (fOpt.sigmaQAccRA)^2], 3));
    Qori = diag([zeros(1,6) repelem((fOpt.sigmaQOriMP)^2, 1, 4) ...
                 zeros(1,6) repelem((fOpt.sigmaQOriLA)^2, 1, 4) ...
                 zeros(1,6) repelem((fOpt.sigmaQOriRA)^2, 1, 4)]);
    Q = G * Q * G' + Qori;
    % initialise covariance in the state estimate
    if islogical(P0) && ~P0
        P_tilde = Q;
    elseif isscalar(P0)
        P_tilde = P0*I_N;
    else
        P_tilde = P0;
    end
    
    nMeasure = 12;
    H = zeros(nMeasure, nStates);
    H(idxMOriMP, idxOriMP) = eye(4, 4);
    H(idxMOriLA, idxOriLA) = eye(4, 4);
    H(idxMOriRA, idxOriRA) = eye(4, 4);

    Rdiag = repelem([(fOpt.sigmaROriMP)^2 (fOpt.sigmaROriLA)^2 (fOpt.sigmaROriRA)^2], 4);
    R = diag(Rdiag);
    
    R_uwb = zeros(3,3);
    R_uwb(1,1) = fOpt.sigmaUwbMPLA.^2;
    R_uwb(2,2) = fOpt.sigmaUwbMPRA.^2;
    R_uwb(3,3) = fOpt.sigmaUwbLARA.^2;
    
    % check that all accelerometer measurements are equal dimensions
    [nSamples, ~] = size(gfrAccMP);
    validateattributes(gfrAccMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(gfrAccRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 3});
    validateattributes(qMP, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qLA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    validateattributes(qRA, {'numeric'}, {'2d', 'nrows', nSamples, 'ncols', 4});
    
    % local variable assignment for readability
    u_k = [gfrAccMP, gfrAccLA, gfrAccRA]';
    y_k = [qMP, qLA, qRA]';
    uwbMPLARA = [uwb_mea.left_tibia_mid_pelvis,...
                 uwb_mea.mid_pelvis_right_tibia,...
                 uwb_mea.left_tibia_right_tibia];
                 
    % allocate memory to store apriori and aposteriori state estimates, xhat,
    % and error covariances in the state estimate, P_pri, P_pos
    xhat_pri = nan(nSamples, nStates);
    P_pri    = nan(nStates, nStates, nSamples);

    xhat_pos = nan(nSamples, nStates);
    P_pos    = nan(nStates, nStates, nSamples);
    
    xhat_con = nan(nSamples, nStates);
    P_con    = nan(nStates, nStates, nSamples);

    debug_dat = struct;
    debug_dat.LFEO = nan(nSamples, 3); debug_dat.RFEO = nan(nSamples, 3);
    debug_dat.LFEP = nan(nSamples, 3); debug_dat.RFEP = nan(nSamples, 3);
    debug_dat.qLTH = nan(nSamples, 4); debug_dat.qRTH = nan(nSamples, 4);

    debug_dat.predState = nan(nSamples, nStates);
    debug_dat.predP = nan(nStates, nStates, nSamples);
%     debug_dat.zuptStateL = false(nSamples, 1);
%     debug_dat.zuptStateR = false(nSamples, 1);
    debug_dat.zuptState = nan(nSamples, nStates);
    debug_dat.zuptP = nan(nStates, nStates, nSamples);
    debug_dat.uwbuptState = nan(nSamples, nStates);
    debug_dat.uwbuptP = nan(nStates, nStates, nSamples);
    debug_dat.cstrState = nan(nSamples, nStates);
    debug_dat.cstrP = nan(nStates, nStates, nSamples);
    debug_dat.cstrStateU = false(nSamples, 1);
    
    debug_dat.zuptStateL = bIsStatLA;
    debug_dat.zuptStateR = bIsStatRA;
    
    refSide = 'N';
    floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
    
    modHundredApplyMeas = mod(fOpt.applyMeas, 100);
    divHundredApplyMeas = idivide(int32(fOpt.applyMeas), 100, 'floor');
    applyCstrW = idivide(int32(fOpt.applyCstr), 100, 'floor');
    
    % applyMeas update orientation initialization
    if fOpt.applyMeas
        idxMVelMP = nMeasure+1:nMeasure+3;
        idxMVelLA = nMeasure+4:nMeasure+6;
        idxMVelRA = nMeasure+7:nMeasure+9;
        nMeasure = nMeasure+9;
        
        H(end+1:end+9, :) = zeros(9, nStates);
        H(idxMVelMP, idxVelMP) = eye(3);
        H(idxMVelLA, idxVelLA) = eye(3);
        H(idxMVelRA, idxVelRA) = eye(3);
        
        Rdiag = diag(R);
        Rdiag(end+1:end+9) = repelem([(fOpt.sigmaZuptMP)^2 ...
            (fOpt.sigmaZuptLA)^2 (fOpt.sigmaZuptRA)^2], 3);
        R = diag(Rdiag);
        
        y_k(end+1:end+9, :) = zeros(9, nSamples);
    end
    
    % applyMeas flat floor assumption initialization
    if modHundredApplyMeas >= 2 % add more zupt features
        floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
        
        switch (modHundredApplyMeas)
            case 2
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 3
                targetL = idxPosLA;
                targetR = idxPosRA;
            case 4
                targetL = idxPosLA;
                targetR = idxPosRA;
            case 5
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 6
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 7
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 8
                targetL = idxPosLA;
                targetR = idxPosRA;
            otherwise % 21,31-33,41-43
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
        end

        targetLN = length(targetL); targetRN = length(targetR);
        targetN = targetLN + targetRN;
        idxMPosLA = nMeasure+1:nMeasure+targetLN;
        idxMPosRA = idxMPosLA(end)+1:idxMPosLA(end)+targetRN;
        nMeasure = nMeasure + targetN;

        H(end+1:end+targetN, :) = zeros(targetN, nStates);
        H(idxMPosLA, targetL) = eye(targetLN);
        H(idxMPosRA, targetR) = eye(targetRN);

        Rdiag = diag(R);
        Rdiag(end+1:end+targetLN) = repelem([fOpt.sigmaRPosLA^2], targetLN);
        Rdiag(end+1:end+targetRN) = repelem([fOpt.sigmaRPosRA^2], targetRN);
        R = diag(Rdiag);

        y_k(end+1:end+targetN, :) = floorZ*ones(targetN, nSamples);  
    end
    
    if (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77)
        % reset covariance of left foot to certain value
        if (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77)
            modTenApplyMeas = 2;
        else
            modTenApplyMeas = mod(fOpt.applyMeas, 10);
        end
        
        switch modTenApplyMeas
            case 1
                target = idxPosLA;
                sigmas = [fOpt.sigmaRPosLARecenter];
                idxMCov1 = nMeasure+1:nMeasure+3;
            case 2
                target = [idxPosMP idxPosLA idxPosRA];
                sigmas = [fOpt.sigmaRPosMPLimit fOpt.sigmaRPosLALimit ...
                          fOpt.sigmaRPosRALimit];
                idxMCov1 = nMeasure+1:nMeasure+9;
            case 3
                target = [idxPosMP idxPosLA idxPosLA idxPosRA];
                sigmas = [fOpt.sigmaRPosMPLimit fOpt.sigmaRPosLALimit ...
                          fOpt.sigmaRPosLARecenter fOpt.sigmaRPosRALimit ];
                idxMCov1 = [nMeasure+1:nMeasure+6 nMeasure+10:nMeasure+12];
                idxMCov2 = [nMeasure+1:nMeasure+3 nMeasure+7:nMeasure+12];
        end
        
        targetN = length(target);
        idx = nMeasure+1:nMeasure+targetN;
        nMeasure = nMeasure+targetN;
        
        H(end+1:end+targetN, :) = zeros(targetN, nStates);
        for i=1:3:targetN
            idx2 = idx(i:i+2);
            target2 = target(i:i+2);
            H(idx2, target2) = eye(3, 3);
        end
        
        y_k(end+1:end+targetN, :) = zeros(targetN, nSamples);
        
        Rdiag = diag(R);
        Rdiag(end+1:end+targetN) = repelem(sigmas.^2, 3);
        R = diag(Rdiag);
    end
    
    if (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77)
        % bit 1: pelvis vel = ankle average vel
        % bit 2: pelvis x y pos = ankle average x y pos
        % bit 3: pelvis z pos = initial pelvis z pos
        modTenApplyMeas = mod(fOpt.applyMeas, 10);
        if bitand(modTenApplyMeas, 1)
            targetN = 3;
            idxMVelMPLARA = nMeasure+1:nMeasure+targetN;
            H(end+1:end+targetN, :) = zeros(targetN, nStates);
            H(idxMVelMPLARA, idxVelMP) = -eye(3,3);
            H(idxMVelMPLARA, idxVelLA) = 0.5*eye(3,3);
            H(idxMVelMPLARA, idxVelRA) = 0.5*eye(3,3);

            Rdiag = diag(R);
            Rdiag(end+1:end+targetN) = repelem(fOpt.sigmaRVelMPLARA.^2, targetN);
            R = diag(Rdiag);

            y_k(end+1:end+targetN, :) = zeros(targetN, nSamples);
            
            nMeasure = nMeasure+targetN;
        end
        if bitand(modTenApplyMeas, 2)
            targetN = 2;
            idxMPosMPLARA = nMeasure+1:nMeasure+targetN;
            H(end+1:end+targetN, :) = zeros(targetN, nStates);
            H(idxMPosMPLARA, idxPosMP(1:2)) = -eye(2,2);
            H(idxMPosMPLARA, idxPosLA(1:2)) = 0.5*eye(2,2);
            H(idxMPosMPLARA, idxPosRA(1:2)) = 0.5*eye(2,2);

            Rdiag = diag(R);
            Rdiag(end+1:end+targetN) = repelem(fOpt.sigmaRPosMPLARA.^2, targetN);
            R = diag(Rdiag);

            y_k(end+1:end+targetN, :) = zeros(targetN, nSamples);
            
            nMeasure = nMeasure+targetN;
        end
        if bitand(modTenApplyMeas, 4)
            idxMPosZMP = nMeasure+1:nMeasure+1;
            nMeasure = nMeasure+1;

            H(end+1:end+1, :) = zeros(1, nStates);
            H(idxMPosZMP, idxPosMP(3)) = 1;

            Rdiag = diag(R);
            Rdiag(end+1:end+1) = fOpt.sigmaRPosZMP;
            R = diag(Rdiag);

            y_k(end+1:end+1, :) = x0(idxPosMP(3));
        end
    end
      
    alphalimit = struct('lkmin', fOpt.alphaLKmin, 'lkmax', fOpt.alphaLKmax, ...
                        'rkmin', fOpt.alphaRKmin, 'rkmax', fOpt.alphaRKmax);
    if (fOpt.applyCstr >= 1 && fOpt.applyCstr <= 8) || ...
        (fOpt.applyCstr >= 71 && fOpt.applyCstr <= 78) || ...
        (fOpt.applyCstr >= 201 && fOpt.applyCstr <= 208) || ...
        (fOpt.applyCstr >= 271 && fOpt.applyCstr <= 278)    
        D = zeros(6, nStates);
        D(1:3,idxPosMP) = -eye(3, 3);
        D(1:3,idxPosLA) = eye(3, 3);
        D(4:6,idxPosMP) = -eye(3, 3);
        D(4:6,idxPosRA) = eye(3, 3);
        
        I_CN = eye(6, 6);
        P_custom = eye(nStates);
        P_custom(idxPosMP,idxPosMP) = 0;
        
        debug_dat.cstrStateU = true(nSamples, 6);
    elseif fOpt.applyCstr >= 11 && fOpt.applyCstr <= 16
        x2cxIdx = [idxPosMP, idxVelMP, idxPosLA, idxVelLA, idxPosRA, idxVelRA];
        nCStates = 18;
        D = zeros(6, nCStates);
        D(1:3,1:3) = -eye(3, 3);
        D(1:3,7:9) = eye(3, 3);
        D(4:6,1:3) = -eye(3, 3);
        D(4:6,13:15) = eye(3, 3);
        
        optimOpt = optimoptions('fmincon', 'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'OptimalityTolerance', fOpt.optimOptimalityTolerance, ...
            'ConstraintTolerance', fOpt.optimConstraintTolerance, ...
            'MaxFunctionEvaluations', fOpt.optimMaxFunctionEvaluations, ...
            'UseParallel', fOpt.optimUseParallel);
    elseif (fOpt.applyCstr >= 21 && fOpt.applyCstr <= 23) || ...
           (fOpt.applyCstr >= 221 && fOpt.applyCstr <= 223)
        D = zeros(6, nStates);
        D(1:3,idxPosLA) = eye(3, 3);
        D(4:6,idxPosRA) = eye(3, 3);
        
        debug_dat.cstrStateU = true(nSamples, 6);
    elseif fOpt.applyCstr >= 31 && fOpt.applyCstr <= 36
        x2cxIdx = [idxPosLA, idxVelLA, idxPosRA, idxVelRA];
        nCStates = 12;
        D = zeros(6, nCStates);
        D(1:3,1:3) = eye(3, 3);
        D(4:6,7:9) = eye(3, 3);
        
        optimOpt = optimoptions('fmincon', 'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'OptimalityTolerance', fOpt.optimOptimalityTolerance, ...
            'ConstraintTolerance', fOpt.optimConstraintTolerance, ...
            'MaxFunctionEvaluations', fOpt.optimMaxFunctionEvaluations, ...
            'UseParallel', fOpt.optimUseParallel);
    elseif fOpt.applyCstr >= 51 && fOpt.applyCstr <= 54
        % 01 constraints + MP/LA/RA zpos = floor zpos
        floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
        
        D = zeros(9, nStates);
        D(1:3,idxPosMP) = -eye(3, 3);
        D(1:3,idxPosLA) = eye(3, 3);
        D(4:6,idxPosMP) = -eye(3, 3);
        D(4:6,idxPosRA) = eye(3, 3);
        D(7:7,idxPosMP(3)) = eye(1, 1);
        D(8:8,idxPosLA(3)) = eye(1, 1);
        D(9:9,idxPosRA(3)) = eye(1, 1);
        
        I_CN = eye(9, 9);
        P_custom = eye(nStates);
        P_custom(idxPosMP,idxPosMP) = 0;
        
        debug_dat.cstrStateU = true(nSamples, 9);
        debug_dat.cstrStateU(:,7:9) = false;
    elseif (fOpt.applyCstr >= 101 && fOpt.applyCstr <= 106) || ...
        (fOpt.applyCstr >= 111 && fOpt.applyCstr <= 116) || ...
        (fOpt.applyCstr >= 121 && fOpt.applyCstr <= 126) || ...
        (fOpt.applyCstr >= 131 && fOpt.applyCstr <= 136)
        optimOpt = optimoptions('fmincon', 'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'OptimalityTolerance', fOpt.optimOptimalityTolerance, ...
            'ConstraintTolerance', fOpt.optimConstraintTolerance, ...
            'MaxFunctionEvaluations', fOpt.optimMaxFunctionEvaluations, ...
            'UseParallel', fOpt.optimUseParallel);
    elseif (fOpt.applyCstr >= 145 && fOpt.applyCstr <= 146) || ...
        (fOpt.applyCstr >= 155 && fOpt.applyCstr <= 156) || ...
        (fOpt.applyCstr >= 165 && fOpt.applyCstr <= 166) || ...
        (fOpt.applyCstr >= 175 && fOpt.applyCstr <= 176)
        P_custom = eye(nStates);
        P_custom(idxPosMP,idxPosMP) = 0;
    end

    for n = 1:nSamples
    %% -----------------------------------------------------------------------
    % Prediction Step using accelerometer measurements ----    
        x_min = F * x_tilde + G * u_k(:,n) ;
        P_min = F * P_tilde * F' + Q;
        xhat_pri(n,:) = x_min;
        P_pri(:,:,n)  = P_min;
        
        debug_dat.predState(n,:) = x_min;
        debug_dat.predP(:,:,n) = P_min;
        
    %% ------------------------------------------------------------------------
    % Measurement update step
    % matrices beginnning with 'H_' are the 'observation matrices' that map
    % the variables in the state estimate vector, xhat, to the measurement
    % domain. In this case we are using
        idx = [idxMOriMP idxMOriLA idxMOriRA];
        if refSide == 'N'
            if bIsStatLA(n), refSide = 'L';
            elseif bIsStatRA(n), refSide = 'R'; end
        elseif ~(bIsStatLA(n) || bIsStatRA(n))
            refSide = 'N';
        end
            
        if modHundredApplyMeas >= 1
            if bIsStatMP(n) idx(end+1:end+3) = idxMVelMP; end
            if bIsStatLA(n) idx(end+1:end+3) = idxMVelLA; end
            if bIsStatRA(n) idx(end+1:end+3) = idxMVelRA; end
        end
        
		if (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77)
            if bIsStatLA(n)
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
            end
            if bIsStatRA(n)
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
            end            
        end
        
        if (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77)
            modTenApplyMeas = mod(fOpt.applyMeas, 10);
            if bitand(modTenApplyMeas, 1)
                idx(end+1:end+3) = idxMVelMPLARA;
            end
            if bitand(modTenApplyMeas, 2)
                idx(end+1:end+2) = idxMPosMPLARA;
            end
            if bitand(modTenApplyMeas, 4)
                idx(end+1:end+1) = idxMPosZMP;
            end
        end
       
        res = y_k(idx, n) - H(idx, :) * x_min;
        K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx));
        x_min1 = x_min + K * res;
        
        if (modHundredApplyMeas == 32) || ...
		   (modHundredApplyMeas >= 70 && modHundredApplyMeas <= 77) || ...
		   (modHundredApplyMeas >= 80 && modHundredApplyMeas <= 87)
            idx2 = [idx idxMCov1];
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        elseif modHundredApplyMeas == 33
            if bIsStatLA(n), idx2 = [idx idxMCov2];
            else, idx2 = [idx idxMCov1]; end
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        else
            idx2 = idx;
        end
        P_min1 = (I_N - K * H(idx2, :)) * P_min;
        
        if fOpt.applyMeas
            debug_dat.zuptState(n,:) = x_min1;
            debug_dat.zuptP(:,:,n) = P_min1;
        end
        
		x_plus = x_min1;
		P_plus = P_min1;
        
        x_plus(idxOriMP,1) = quatnormalize(x_plus(idxOriMP,1)')';
        x_plus(idxOriLA,1) = quatnormalize(x_plus(idxOriLA,1)')';
        x_plus(idxOriRA,1) = quatnormalize(x_plus(idxOriRA,1)')';
        
        xhat_pos(n, :) = x_plus;
        P_pos(:, :, n)  = P_plus;
        
    %% -----------------------------------------------------------------------
    % Constraint update step ---- 
        PELV_CS = quat2rotm(x_plus(idxOriMP,1)');
        LTIB_CS = quat2rotm(x_plus(idxOriLA,1)');
        RTIB_CS = quat2rotm(x_plus(idxOriRA,1)');
        % Test frankenstein constraint
        if (applyCstrW >= 3 && applyCstrW <= 5)
            applyCstrX = mod(idivide(int32(fOpt.applyCstr), 10, 'floor'), 10);
            applyCstrModTen = mod(fOpt.applyCstr, 10);
            
            sckfAlpha = fOpt.sckfAlpha;
            sckfThreshold = fOpt.sckfThreshold;
            if applyCstrW == 3 || applyCstrW == 5
                x_tilde = x_plus;
                P_tilde = P_plus;
                
                idxCPosMP = idxPosMP;
                idxCPosLA = idxPosLA;
                idxCPosRA = idxPosRA;
                I_N2 = I_N;
                
                if applyCstrW == 3
                    DbaseRowN = 4;
                else
                    DbaseRowN = 7;
                end
            else
                % applyCstrW == 4
                idx = [idxPosMP idxPosLA idxPosRA];
                x_tilde = x_plus(idx);
                P_tilde = P_plus(idx, idx);
                
                idxCPosMP = 1:3;
                idxCPosLA = 4:6;
                idxCPosRA = 7:9;
                I_N2 = eye(length(x_tilde), length(x_tilde));
                
                DbaseRowN = 4;
            end
            
            % preprocessing for knee angle inequality
            % additional process to ensure knee angle does not increase
            % during iteration
            if (applyCstrModTen >= 5 && applyCstrModTen <= 6)
                LFEM_z = x_tilde(idxCPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxCPosLA,1);
                RFEM_z = x_tilde(idxCPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxCPosRA,1);
                     
                alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                                 -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
                alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                                 -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
                alphalimit.lkmax = min(alpha_lk, fOpt.alphaLKmax);
                alphalimit.rkmax = min(alpha_rk, fOpt.alphaRKmax);
            end
            
            tmpLKmin = sin(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,1) - ...
                       cos(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,3);
            tmpLKmax = sin(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,1) - ...
                       cos(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,3);
            tmpRKmin = sin(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,1) - ...
                       cos(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,3);
            tmpRKmax = sin(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,1) - ...
                       cos(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,3);
                   
            for i=0:fOpt.sckfMaxIter
                % Step 1: construct D matrix              
                % calculate the z axis of femur and tibia
                LFEM_z = x_tilde(idxCPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxCPosLA,1);
                RFEM_z = x_tilde(idxCPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxCPosRA,1);
                g_dlfem = norm(LFEM_z, 2);
                g_drfem = norm(RFEM_z, 2);
                
                % D matrix construction: femur length constraint and knee hinge joint
                D = zeros(DbaseRowN, length(x_tilde));
                res = zeros(DbaseRowN, 1);
                sigmaCstr = zeros(DbaseRowN, 1);
                D(1, idxCPosMP) = LFEM_z'/g_dlfem;
                D(1, idxCPosLA) = -LFEM_z'/g_dlfem;
                D(2, idxCPosMP) = RFEM_z'/g_drfem;
                D(2, idxCPosRA) = -RFEM_z'/g_drfem;
                D(3, idxCPosMP) = LTIB_CS(:,2)';
                D(3, idxCPosLA) = -LTIB_CS(:,2)';
                D(4, idxCPosMP) = RTIB_CS(:,2)';
                D(4, idxCPosRA) = -RTIB_CS(:,2)';
                
                res(1, 1) = dLFemur - g_dlfem;
                res(2, 1) = dRFemur - g_drfem;
                res(3, 1) = (-dPelvis/2*PELV_CS(:,2) + ...
                    dLTibia*LTIB_CS(:,3))'*LTIB_CS(:,2) - D(3,:)*x_tilde;
                res(4, 1) = (dPelvis/2*PELV_CS(:,2) + ...
                    dRTibia*RTIB_CS(:,3))'*RTIB_CS(:,2) - D(4,:)*x_tilde;
                        
                if applyCstrW == 5 % apply uwb at cstr update
                    diffMPLA = x_tilde(idxPosMP)'-x_tilde(idxPosLA)';
                    diffMPRA = x_tilde(idxPosMP)'-x_tilde(idxPosRA)';
                    diffLARA = x_tilde(idxPosLA)'-x_tilde(idxPosRA)';
                    % the observation model
                    hUwbEst = [vecnorm(diffMPLA, 2, 2);
                                 vecnorm(diffMPRA, 2, 2);
                                 vecnorm(diffLARA, 2, 2)];
            
                    D(5, idxCPosMP) = diffMPLA./hUwbEst(1);
                    D(5, idxCPosLA) = -diffMPLA./hUwbEst(1);
                    D(6, idxCPosMP) = diffMPRA./hUwbEst(2);
                    D(6, idxCPosRA) = -diffMPRA./hUwbEst(2);
                    D(7, idxCPosLA) = diffLARA./hUwbEst(3);
                    D(7, idxCPosRA) = -diffLARA./hUwbEst(3);

                    res(5:7, 1) = uwbMPLARA(n,:)' - hUwbEst;
                    sigmaCstr(5, 1) = fOpt.sigmaUwbMPLA;
                    sigmaCstr(6, 1) = fOpt.sigmaUwbMPRA;
                    sigmaCstr(7, 1) = fOpt.sigmaUwbLARA;
                end
                

                  
                % add knee inequality constraint via active set
                alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                                 -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
                alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                                 -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
                if alpha_lk < alphalimit.lkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpLKmin';
                    D(end, idxCPosLA) = -tmpLKmin';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmin ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                elseif alpha_lk > alphalimit.lkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpLKmax';
                    D(end, idxCPosLA) = -tmpLKmax';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmax ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                end
                    
                if alpha_rk < alphalimit.rkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpRKmin';
                    D(end, idxCPosRA) = -tmpRKmin';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmin ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                elseif alpha_rk > alphalimit.rkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxCPosMP) = tmpRKmax';
                    D(end, idxCPosRA) = -tmpRKmax';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmax ...
                        - D(end, :)*x_tilde;
                    sigmaCstr(end+1, 1) = 0;
                end
                
                % add flat floor constraint
                if (applyCstrX >= 1 && applyCstrX <= 3)
                    if refSide == 'L'
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxCPosLA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        sigmaCstr(end+1, 1) = 0;
                    elseif refSide == 'R'
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxCPosRA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        sigmaCstr(end+1, 1) = 0;
                    end
                elseif (applyCstrX >= 7 && applyCstrX <= 9)
%                     stepCstrCount = 0;
%                     if bIsStatLA(n) || x_tilde(idxCPosLA(3), 1) < floorZ
                    if refSide == 'L' || x_tilde(idxCPosLA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxCPosLA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        sigmaCstr(end+1, 1) = 0;
%                         stepCstrCount = stepCstrCount + 1;
%                     end
%                     if bIsStatRA(n) || x_tilde(idxCPosRA(3), 1) < floorZ
                    elseif refSide == 'R' || x_tilde(idxCPosRA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxCPosRA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        sigmaCstr(end+1, 1) = 0;
%                         stepCstrCount = stepCstrCount + 1;
                    end
%                     if stepCstrCount == 2
%                         P_tilde_cstr = D*P_tilde*D';
%                         % if condition number is high enough (singular
%                         % matrix), delete last constraint
%                         if cond(P_tilde_cstr) > 1e5
%                             D(end, :) = [];
%                             res(end, :) = [];
%                         end
%                     end
                end
                
                % add static foot pos x y on foot step constraint
                if (applyCstrX == 2) || (applyCstrX == 8)
%                     if bIsStatLA(n)
                    if refSide == 'L'
                        D(end+1:end+2, :) = zeros(2, length(x_tilde));
                        D(end-1:end, idxCPosLA(1:2)) = eye(2,2);
                        res(end+1:end+2, 1) = x_plus(idxCPosLA(1:2), 1) - x_tilde(idxCPosLA(1:2), 1);
                        sigmaCstr(end+1:end+2, 1) = 0;
%                     elseif bIsStatRA(n)
                    elseif refSide == 'R'
                        D(end+1:end+2, :) = zeros(2, length(x_tilde));
                        D(end-1:end, idxCPosRA(1:2)) = eye(2,2);
                        res(end+1:end+2, 1) = x_plus(idxCPosRA(1:2), 1) - x_tilde(idxCPosRA(1:2), 1);
                        sigmaCstr(end+1:end+2, 1) = 0;
                    end
                end
                
                if (applyCstrModTen >= 3 && applyCstrModTen <= 4)
                    P_custom = eye(length(x_tilde), length(x_tilde));
                else
                    P_custom = P_tilde;
                end
                
                if (applyCstrX == 3) || (applyCstrX == 9)
%                     if bIsStatLA(n)
                    if refSide == 'L'
                        P_custom(idxCPosLA(1:2), idxCPosLA(1:2)) = 0;
%                     elseif bIsStatRA(n)
                    elseif refSide == 'R'
                        P_custom(idxCPosRA(1:2), idxCPosRA(1:2)) = 0;
                    end
                end
                
                % Step 3: calculate Ri               
                Ri = sckfAlpha*D*P_tilde*D'*exp(-i);
                Si_denom = diag(D * P_tilde* D');
                Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ Si_denom;
                Si_ignoreindex = (Si_denom < sigmaCstr) | isnan(Si) | (Si_denom < 1e-10);
                Si(Si_ignoreindex) = sckfThreshold+1;
                if sum(Si < sckfThreshold) == 0
                    break
                end

                % Step 4: update x_tilde and P_tilde
                switch mod(fOpt.applyCstr, 10)
                    case 1
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                        P_tilde = (I_N2-Kk*D)*P_tilde*(I_N2-Kk*D)' + Kk*Ri*Kk';
                    case 2
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                    case 3
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
%                         Kk = D'*(D*D'+Ri)^(-1);
                        P_tilde = (I_N2-Kk*D)*P_tilde*(I_N2-Kk*D)' + Kk*Ri*Kk';
                    case 4
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
%                         Kk = D'*(D*D'+Ri)^(-1);
                    case 5
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                        P_tilde = (I_N2-Kk*D)*P_tilde*(I_N2-Kk*D)' + Kk*Ri*Kk';
                    case 6
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                    otherwise
                        Kk = 0;
                end
                
                dx = Kk*(res);
                x_tilde = x_tilde + dx;
            end
			
            if applyCstrW == 3 || applyCstrW == 5
                P_tilde = P_plus;
            else % applyCstrW == 4
                x_tilde2 = x_tilde;
                x_tilde = x_plus;
                x_tilde(idx) = x_tilde2;
                P_tilde = P_plus;
            end
        
		end
            
            x0 = x_plus(idx);
            
            x = fmincon(@(x) L2(x, x0, W, 20), ...
                x0, [], [], [], [], [], [], nonlcon, optimOpt);
            
            % normalize the quaternions
            x_tilde = x_plus;
            x_tilde(idx) = x;
            x_tilde(idxOriMP) = quatnormalize(x_tilde(idxOriMP, 1)');
            x_tilde(idxOriLA) = quatnormalize(x_tilde(idxOriLA, 1)');
            x_tilde(idxOriRA) = quatnormalize(x_tilde(idxOriRA, 1)');
            P_tilde = P_plus;
        end
        
        xhat_con(n, :) = x_tilde;
        P_con(:, :, n)  = P_tilde;
        debug_dat.cstrState(n,:) = x_tilde;
        debug_dat.cstrP(:,:,n) = P_tilde;
        
        debug_dat.LFEO(n, :) = x_tilde(idxPosLA) + dLTibia * LTIB_CS(:, 3);
        debug_dat.RFEO(n, :) = x_tilde(idxPosRA) + dRTibia * RTIB_CS(:, 3);
        debug_dat.LFEP(n, :) = x_tilde(idxPosMP) + dPelvis/2 * PELV_CS(:, 2);
        debug_dat.RFEP(n, :) = x_tilde(idxPosMP) - dPelvis/2 * PELV_CS(:, 2);
        
%         if fOpt.applyCstr
        LFEM_z = (debug_dat.LFEP(n,:)-debug_dat.LFEO(n,:))'; 
        LFEM_y = LTIB_CS(:,2);
        LFEM_x = cross(LFEM_y, LFEM_z);
        LFEM_z = LFEM_z / norm(LFEM_z);
        LFEM_y = LFEM_y / norm(LFEM_y);
        LFEM_x = LFEM_x / norm(LFEM_x);
        RFEM_z = (debug_dat.RFEP(n,:)-debug_dat.RFEO(n,:))';
        RFEM_y = RTIB_CS(:,2);
        RFEM_x = cross(RFEM_y, RFEM_z);
        RFEM_y = RFEM_y / norm(RFEM_y);
        RFEM_z = RFEM_z / norm(RFEM_z);
        RFEM_x = RFEM_x / norm(RFEM_x);
        debug_dat.qLTH(n, :) = rotm2quat([LFEM_x LFEM_y LFEM_z]);
        debug_dat.qRTH(n, :) = rotm2quat([RFEM_x RFEM_y RFEM_z]);
%         end
    end
    debug_dat.y_k = y_k';
end

function d = solveKneeAngleFromDist(PELV_CS, LTIB_CS, RTIB_CS, ...
                            dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, ...
                            dLUwbDist, dRUwbDist)
    lbuf = (dLTibia^2+dLFemur^2-dLUwbDist^2)/(2*dLTibia*dLFemur);
    rbuf = (dRTibia^2+dRFemur^2-dRUwbDist^2)/(2*dRTibia*dRFemur);
    lbuf = min(max(lbuf, -1), 1);
    rbuf = min(max(rbuf, -1), 1);
    
    % calculate alpha_lk and alpha_rk
    alpha_lk = pi - acos(lbuf);
    alpha_rk = pi - acos(rbuf);

    % setup the constraint equations
    d = [ (dPelvis/2*PELV_CS(:,2) ...
             -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
             +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
             -dLTibia*LTIB_CS(:,3)) ; ...
            (-dPelvis/2*PELV_CS(:,2)+ ...
             -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
             +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
             -dRTibia*RTIB_CS(:,3)) ];
end

function d = solve_linhjc_d(pMP, pLA, pRA, PELV_CS, LTIB_CS, RTIB_CS, ...
                            dPelvis, dLFemur, dRFemur, dLTibia, dRTibia)
    LKNE = pLA + dLTibia*LTIB_CS(:,3);
    RKNE = pRA + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = pMP+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = pMP-dPelvis/2*PELV_CS(:,2)-RKNE;

    % calculate the z axis of the tibia
%     LTIB_z = LTIB_CS(:,3);
%     RTIB_z = RTIB_CS(:,3);
    LFEM_z__TIB = LTIB_CS\LFEM_z;
    RFEM_z__TIB = RTIB_CS\RFEM_z;

    % calculate alpha_lk and alpha_rk
    alpha_lk = atan2(-LFEM_z__TIB(3), -LFEM_z__TIB(1)) + 0.5*pi;
    alpha_rk = atan2(-RFEM_z__TIB(3), -RFEM_z__TIB(1)) + 0.5*pi;
%     alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
%     alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));

    % setup the constraint equations
    d = [ (dPelvis/2*PELV_CS(:,2) ...
             -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
             +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
             -dLTibia*LTIB_CS(:,3)) ; ...
            (-dPelvis/2*PELV_CS(:,2)+ ...
             -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
             +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
             -dRTibia*RTIB_CS(:,3)) ];
end

function d = solve_linhjc_kac_d(pMP, pLA, pRA, PELV_CS, LTIB_CS, RTIB_CS, ...
                            dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, ...
                            alphalimit)
    LKNE = pLA + dLTibia*LTIB_CS(:,3);
    RKNE = pRA + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = pMP+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = pMP-dPelvis/2*PELV_CS(:,2)-RKNE;

    % normalize z-axis of femur
    LFEM_z = LFEM_z/norm(LFEM_z);
    RFEM_z = RFEM_z/norm(RFEM_z);

    % _TIB_CS is _TIB_CS described in world frame, or rotm from tib2world frame
    % therefore, inverse is from wrold to tib frame
    LFEM_z__TIB = LTIB_CS\LFEM_z;
    RFEM_z__TIB = RTIB_CS\RFEM_z;

    %global alpha_lk alpha_rk 
    alpha_lk = atan2(-LFEM_z__TIB(3), -LFEM_z__TIB(1)) + 0.5*pi;
    alpha_rk = atan2(-RFEM_z__TIB(3), -RFEM_z__TIB(1)) + 0.5*pi;
    
    alpha_lk = max(alphalimit.lkmin, min(alpha_lk, alphalimit.lkmax));
    alpha_rk = max(alphalimit.rkmin, min(alpha_rk, alphalimit.rkmax));
    
    % setup the constraint equations
    d = [ (dPelvis/2*PELV_CS(:,2) ...
             -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
             +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
             -dLTibia*LTIB_CS(:,3)) ; ...
            (-dPelvis/2*PELV_CS(:,2)+ ...
             -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
             +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
             -dRTibia*RTIB_CS(:,3)) ];
end

function [la ra] = calcKneeAngle(x, PELV_CS, LTIB_CS, RTIB_CS, ...
                        dPelvis, dLTibia, dRTibia)
    idxPosMP = 1:3;
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation

    LKNE = x(idxPosLA) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA) + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP)-dPelvis/2*PELV_CS(:,2)-RKNE;

    % normalize z-axis of femur
    LFEM_z = LFEM_z/norm(LFEM_z);
    RFEM_z = RFEM_z/norm(RFEM_z);

    % _TIB_CS is _TIB_CS described in world frame, or rotm from tib2world frame
    % therefore, inverse is from wrold to tib frame
    LFEM_z__TIB = LTIB_CS\LFEM_z;
    RFEM_z__TIB = RTIB_CS\RFEM_z;

    %global alpha_lk alpha_rk 
    la = atan2(-LFEM_z__TIB(3), -LFEM_z__TIB(1)) + 0.5*pi;
    ra = atan2(-RFEM_z__TIB(3), -RFEM_z__TIB(1)) + 0.5*pi;
end

function y = L2Dist(x, x0, S, exponent)
%x^2 is monotomically increasing at any point not at 0, so traditional
%L2 norm involving sqrt is unnecessary, can use x^2 to find same
%location of min cost in constrained region with less computational
%cost.
%using inverse of covariance rather than I will make state est over time
%less smooth but more accurate over the average of the interval
    qW = norm(x0); %weighting for q deviation to account for diffin units between pos and q.
    res = (x-x0);
    qIdx = [7:10 17:20 27:30]';
    res(qIdx) = qW*res(qIdx);
    res = S\res; %add res*inv(S) to scale cost by certainty
    %have also tried 2,4,6,8,14,100,1000 around >= 14 greatly increases speed of finding solution, not much difference etween 100 and 1000
    y = (res'*res)^exponent;
end

function y = L2(x, x0, S, exponent)
    res = (x-x0);
    res = S\res; %add res*inv(S) to scale cost by certainty
    y = (res'*res)^exponent;
end

%% constraint function for cstr 101 - 106
function [c, ceq] = hjc_nonlcon(x, dPelvis, dLFemur, dRFemur, dLTibia, dRTibia)
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    % idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    % idxVelLA = 14:16; % column idx corresponding to the left ankle velocity
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    % idxVelRA = 24:26; % column idx corresponding to the right ankle velocity
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation

    LTIB_CS = quat2rotm(x(idxOriLA,1)');
    RTIB_CS = quat2rotm(x(idxOriRA,1)');
    PELV_CS = quat2rotm(x(idxOriMP,1)');

    LKNE = x(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA,1) + dRTibia*RTIB_CS(:,3);
    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;
    
    % normalize z-axis of femur
    LFEM_z = LFEM_z/norm(LFEM_z);
    RFEM_z = RFEM_z/norm(RFEM_z);
    
    % calculate the z axis of the tibia
    LTIB_z = LTIB_CS(:,3);
    RTIB_z = RTIB_CS(:,3);

    % calculate alpha_lk and alpha_rk
    alpha_lk = acos(dot(LFEM_z, LTIB_z)/(norm(LFEM_z)*norm(LTIB_z)));
    alpha_rk = acos(dot(RFEM_z, RTIB_z)/(norm(RFEM_z)*norm(RTIB_z)));

    % setup the constraint equations
    ceq = [x(idxPosLA) - x(idxPosMP) - (dPelvis/2*PELV_CS(:,2) ...
        -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
        +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
        -dLTibia*LTIB_CS(:,3)) ; ...
        x(idxPosRA) - x(idxPosMP) - (-dPelvis/2*PELV_CS(:,2)+ ...
        -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
        +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
        -dRTibia*RTIB_CS(:,3)) ];

    c = [];
end

%% constraint function for cstr 111 - 116
function [c, ceq] = hjc_nonlcon_kac(x, dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, alphalimit)
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    % idxVelMP = 4:6; % column idx corresponding to the mid-pelvis velocity
	idxOriMP = 7:10; % column idx corresponding to the mid-pelvis orientation
    idxPosLA = 11:13; % column idx corresponding to the left ankle position
    % idxVelLA = 14:16; % column idx corresponding to the left ankle velocity
    idxOriLA = 17:20; % column idx corresponding to the left ankle orientation
    idxPosRA = 21:23; % column idx corresponding to the right ankle position
    % idxVelRA = 24:26; % column idx corresponding to the right ankle velocity
    idxOriRA = 27:30; % column idx corresponding to the right ankle orientation

    LTIB_CS = quat2rotm(x(idxOriLA,1)');
    RTIB_CS = quat2rotm(x(idxOriRA,1)');
    PELV_CS = quat2rotm(x(idxOriMP,1)');

    LKNE = x(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA,1) + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;

    % normalize z-axis of femur
    LFEM_z = LFEM_z/norm(LFEM_z);
    RFEM_z = RFEM_z/norm(RFEM_z);
    
    % _TIB_CS is _TIB_CS described in world frame, or rotm from tib2world frame
    % therefore, inverse is from wrold to tib frame
    LFEM_z__TIB = LTIB_CS\LFEM_z;
    RFEM_z__TIB = RTIB_CS\RFEM_z;

    %global alpha_lk alpha_rk 
    alpha_lk = atan2(-LFEM_z__TIB(3), -LFEM_z__TIB(1)) + 0.5*pi;
    alpha_rk = atan2(-RFEM_z__TIB(3), -RFEM_z__TIB(1)) + 0.5*pi;

%     if alpha_lk < 0, alpha_lk = 0;
%     elseif alpha_lk > pi, alpha_lk = pi; end
%     if alpha_rk < 0, alpha_rk = 0;
%     elseif alpha_rk > pi, alpha_rk = pi; end
    
    % setup the constraint equations
    ceq = [x(idxPosLA) - x(idxPosMP) - (dPelvis/2*PELV_CS(:,2) ...
        -dLFemur*cos(alpha_lk)*LTIB_CS(:,3) ...
        +dLFemur*sin(alpha_lk)*LTIB_CS(:,1) ...
        -dLTibia*LTIB_CS(:,3)) ; ...
        x(idxPosRA) - x(idxPosMP) - (-dPelvis/2*PELV_CS(:,2)+ ...
        -dRFemur*cos(alpha_rk)*RTIB_CS(:,3) ...
        +dRFemur*sin(alpha_rk)*RTIB_CS(:,1) ...
        -dRTibia*RTIB_CS(:,3)) ];

    c = [alphalimit.lkmin - alpha_lk; 
         alpha_lk - alphalimit.lkmax;
         alphalimit.rkmin - alpha_rk; 
         alpha_rk - alphalimit.rkmax];
end

%% constraint function for cstr 121 - 126
function [c, ceq] = hjc_nonlcon_poska(x, PELV_CS, LTIB_CS, RTIB_CS, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia)
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    idxPosLA = 4:6; % column idx corresponding to the left ankle position
    idxPosRA = 7:9; % column idx corresponding to the right ankle position
    
    LKNE = x(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA,1) + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;
    
    % setup the constraint equations
    ceq = [norm(LFEM_z, 2) - dLFemur;
        norm(RFEM_z, 2) - dRFemur;
        dot(LFEM_z, LTIB_CS(:,2)); 
        dot(RFEM_z, RTIB_CS(:,2)) ];

    c = [];
end

%% constraint function for cstr 131 - 136
function [c, ceq] = hjc_nonlcon_kac_poska(x, PELV_CS, LTIB_CS, RTIB_CS, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, alphalimit)
    idxPosMP = 1:3; % column idx corresponding to the mid-pelvis position
    idxPosLA = 4:6; % column idx corresponding to the left ankle position
    idxPosRA = 7:9; % column idx corresponding to the right ankle position
    
    LKNE = x(idxPosLA,1) + dLTibia*LTIB_CS(:,3);
    RKNE = x(idxPosRA,1) + dRTibia*RTIB_CS(:,3);

    % calculate the z axis of the femur
    LFEM_z = x(idxPosMP,1)+dPelvis/2*PELV_CS(:,2)-LKNE;
    RFEM_z = x(idxPosMP,1)-dPelvis/2*PELV_CS(:,2)-RKNE;
    
    % _TIB_CS is _TIB_CS described in world frame, or rotm from tib2world frame
    % therefore, inverse is from wrold to tib frame
    LFEM_z__TIB = LTIB_CS\(LFEM_z/norm(LFEM_z));
    RFEM_z__TIB = RTIB_CS\(RFEM_z/norm(RFEM_z));

    %global alpha_lk alpha_rk 
    alpha_lk = atan2(-LFEM_z__TIB(3), -LFEM_z__TIB(1)) + 0.5*pi;
    alpha_rk = atan2(-RFEM_z__TIB(3), -RFEM_z__TIB(1)) + 0.5*pi;
    
    % setup the constraint equations
    ceq = [norm(LFEM_z) - dLFemur;
        norm(RFEM_z) - dRFemur;
        dot(LFEM_z, LTIB_CS(:,2)); 
        dot(RFEM_z, RTIB_CS(:,2)) ];

    c = [alphalimit.lkmin - alpha_lk; 
         alpha_lk - alphalimit.lkmax;
         alphalimit.rkmin - alpha_rk; 
         alpha_rk - alphalimit.rkmax];
end