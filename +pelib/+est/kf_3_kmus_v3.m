function [ xhat_pri, xhat_con, debug_dat ] = kf_3_kmus_v3(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, ...
    gfrAccLA, bIsStatLA, qLA, ...
    gfrAccRA, bIsStatRA, qRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)
% KF_3_KMUS Kalman Filter for performing sensor fusion on the trajectory of
% three KMUs presumably worn on the body in the following configuration: mid
% pelvis, left ankle, right ankle
% In this state space model, the position and velocity of each kinematic
% measurement unit (KMU) is estimated in 3D space by combining the
% information from each KMU in a kalman filter. NOTE: pay special attention 
% to units:
% position (meters)
% velocity (m/s)
% acceleration (m/2^2)
% uwb_mea (meters)
%
% Author: Luke Wicent Sy, Michael Del Rosario
%
% Inputs::
%   fs - sampling frequency of the magnetic and inertial measurement units
%   sigma_acc - user specified process noise, i.e., the standard deviation
%               in the accelerometer measurements when subjected to a known
%               acceleration
%   x0        - the initial state in the GFR
%   gfrAccMP - the acceleration of the mid-pelvis in the GFR
%   gfrAccLA - the acceleration of the left ankle in the GFR
%   gfrAccRA - the acceleration of the right ankle in the GFR
%   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%                performed by using psuedo-zero velocity measurements 
%   qMP       - mid  pelvis orientation in the GFR (quaternion)
%   qLA       - left  ankle orientation in the GFR (quaternion)
%   qRA       - right ankle orientation in the GFR (quaternion)
%   dPelvis   - pelvis width
%   dRFemur   - right femur length
%   dLFemur   - left femur length
%   dRTibia   - right tibia length
%   dLTibia   - left tibia length
%   uwb_mea    - a structure containing the range measurements (m) between
%   options   - struct containing the ff. settings:
%       applyMeas - turn on/off zero velocity update. boolean
%           001: standard zupt
%           002: standard zupt + floor assumption
%           003: standard zupt + floor assumption + reset at both foot
%           004: standard zupt + floor assumption + series reset at both foot
%           005: standard zupt + floor assumption + zero first floor foot
%           006: standard zupt + floor assumption + pelvis z pos assumption
%           007: standard zupt + floor assumption + pelvis z pos assumption + 
%                zero first floor foot
%           008: standard zupt + floor assumption + reset first floor foot
%           021: standard zupt + floor assumption + estimate projection (W=P^-1) assuming perfect orientation
%           031: standard zupt + floor assumption + recenter at left foot step
%                + one big cov update
%           032: standard zupt + floor assumption + limit pos max covariance
%                + one big cov update
%           033: standard zupt + floor assumption + recenter at left foot step
%                + limit pos max covariance + one big cov update
%           041: standard zupt + floor assumption + recenter at left foot step
%                + series cov update
%           042: standard zupt + floor assumption + limit pos max covariance
%                + series cov update
%           043: standard zupt + floor assumption + recenter at left foot step
%                + limit pos max covariance + series cov update
%           051: standard zupt + floor assumption + pelvis = ankle speed
%           052: standard zupt + floor assumption + pelvis = ankle XY pos
%           053: standard zupt + floor assumption + pelvis = ankle speed
%                + pelvis = ankle XY pos
%           101: 3 point dist + standard zupt
%           102: 3 point dist + standard zupt + floor assumption
%           103: 3 point dist + standard zupt + floor assumption + reset at both foot
%           104: 3 point dist + standard zupt + floor assumption + series reset at both foot
%           105: 3 point dist + standard zupt + floor assumption + zero first floor foot
%           106: 3 point dist + standard zupt + floor assumption + pelvis z pos assumption
%           107: 3 point dist + standard zupt + floor assumption + pelvis z pos assumption + 
%                zero first floor foot
%           121: 3 point dist + standard zupt + floor assumption + estimate projection (W=P^-1) assuming perfect orientation
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
%           351: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update
%           352: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update
%           353: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update
%           354: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update
%           355: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%           356: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P updates
%           371: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor
%           372: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor
%           373: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor
%           374: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor
%           375: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%                + foot step (z pos) = floor
%           376: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P updates
%                + foot step (z pos) = floor
%           381: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor + static foot step (x,y pos)
%           382: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor + static foot step (x,y pos)
%           383: smoothly constraint kf (W=I, early stop)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor + static foot step (x,y pos)
%           384: smoothly constraint kf (W=I, use maxIter)
%                + adj knee + knee ineq + no P update
%                + foot step (z pos) = floor + static foot step (x,y pos)
%           385: smoothly constraint kf (W=P^-1, early stop)
%                + adj knee that only increases + knee ineq + no P update
%                + foot step (z pos) = floor + static foot step (x,y pos)
%           386: smoothly constraint kf (W=P^-1, use maxIter)
%                + adj knee that only increases + knee ineq + no P updates
%                + foot step (z pos) = floor + static foot step (x,y pos)
%         'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
    fOpt = struct('fs', 60, 'applyMeas', false, 'applyUwb', false, ...
        'applyAccBias', false, 'applyCstr', 0, ...
        'sigmaQAccMP', 0.5, 'sigmaQAccLA', 0.5, 'sigmaQAccRA', 0.5, ...
        'sigmaQOriMP', 1e3, 'sigmaQOriLA', 1e3, 'sigmaQOriRA', 1e3, ...
        'sigmaROriMP', 1e-3, 'sigmaROriLA', 1e-3, 'sigmaROriRA', 1e-3, ...
        'sigmaRPosLA', 1e-4, 'sigmaRPosRA', 1e-4, 'sigmaRPosMP', 1e-4, ...
        'sigmaRPosMPLimit', 1e2, 'sigmaRPosLALimit', 1e2, ...
        'sigmaRPosRALimit', 1e2, 'sigmaRPosLARecenter', 1e-3, ...
        'sigmaRVelMPLARA', 1e1, 'sigmaRPosMPLARA', 1e1, ...
        'sigmaCPos', 1e-2, ...
        'sigmaUwbMPLA', 0.2, 'sigmaUwbMPRA', 0.2, 'sigmaUwbLARA', 0.1, ...
        'sigmaZuptMP', 1e-1, 'sigmaZuptLA', 1e-1, 'sigmaZuptRA', 1e-1, ...
        'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
        'alphaRKmin', 0, 'alphaRKmax', pi*8/9, ...
        'optimOptimalityTolerance', 1e-2, ...
        'optimConstraintTolerance', 1e-2, ...
        'optimMaxFunctionEvaluations', 1500, 'optimUseParallel', false, ...
        'sckfAlpha', 0.1, 'sckfThreshold', 100, 'sckfMaxIter', 100);
    
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
    
    % applyMeas set pelvis z position bias to certain height initialization
    if mod(fOpt.applyMeas, 100) >= 6 && mod(fOpt.applyMeas, 100) <= 7
        idxMPosMP = nMeasure+1:nMeasure+1;
        nMeasure = nMeasure+1;
        
        H(end+1:end+1, :) = zeros(1, nStates);
        H(idxMPosMP, idxPosMP(3)) = 1;
        
        Rdiag = diag(R);
        Rdiag(end+1:end+1) = fOpt.sigmaRPosMP;
        R = diag(Rdiag);
        
        y_k(end+1:end+1, :) = x0(idxPosMP(3));
    end
    
    % applyMeas flat floor assumption initialization
    if mod(fOpt.applyMeas, 100) >= 2 % add more zupt features
        floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
        
        switch (mod(fOpt.applyMeas, 100))
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
                refSide = 'N';
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 6
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 7
                refSide = 'N';
                targetL = idxPosLA(3);
                targetR = idxPosRA(3);
            case 8
                refSide = 'N';
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
    
    % applyMeas special measurement update
    if mod(fOpt.applyMeas, 100) == 21
        idxMLHjc = nMeasure+1:nMeasure+3;
        idxMRHjc = nMeasure+4:nMeasure+6;
        nMeasure = nMeasure + 6;
        
        H(end+1:end+6, :) = zeros(6, nStates);
        H(idxMLHjc, idxPosMP) = -eye(3, 3);
        H(idxMLHjc, idxPosLA) = eye(3, 3);
        H(idxMRHjc, idxPosMP) = -eye(3, 3);
        H(idxMRHjc, idxPosRA) = eye(3, 3);
        
        Rdiag = diag(R);
        Rdiag(end+1:end+6) = repelem((fOpt.sigmaCPos)^2, 6);
        R = diag(Rdiag);
        
        y_k(end+1:end+6, :) = zeros(6, nSamples);
    end
    if (mod(fOpt.applyMeas, 100) >= 31 && mod(fOpt.applyMeas, 100) <= 33) || ...
       (mod(fOpt.applyMeas, 100) >= 41 && mod(fOpt.applyMeas, 100) <= 43)
        % reset covariance of left foot to certain value       
        switch (mod(fOpt.applyMeas, 10))
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
        Rdiag = diag(R);
        Rdiag(end+1:end+targetN) = repelem(sigmas.^2, 3);
        R = diag(Rdiag);
    end
    
    if (mod(fOpt.applyMeas, 100) >= 51 && mod(fOpt.applyMeas, 100) <= 53)
        % pelvis average speed / pos = ankle average speed / pos
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
    end
    
    if fOpt.applyMeas >= 100 && fOpt.applyMeas < 200
        idxM3Dist = nMeasure+1:nMeasure+3;
        nMeasure = nMeasure + 3;

        % pseudo UWB measurements corresponding to the euclidean distance between
        % pairs of KMUs.
        % NOTE: the column order of these measurements is important. Assume the
        % following order unless stated otherwise:
        %  uwb_MP_LA_RA =  ['mid-pelvis to left ankle',
        %                   'mid-pelvis to right ankle',
        %                   'left ankle to right ankle'];
        uwb_MP_LA_RA = [uwb_mea.left_tibia_mid_pelvis,...
                        uwb_mea.mid_pelvis_right_tibia,...
                        uwb_mea.left_tibia_right_tibia];

        H(end+1:end+3, :) = zeros(3, nStates);
        
        y_k(end+1:end+3, :) = zeros(3, nSamples);
        % specify the measurement noise in the UWB measurements, these may be
        % different for each KMU sensor pair. It is likely that the range between
        % feet/ankles will be the most accurate due to less "no line of sight"
        % periods. Note: units on sigma_uwb = meters
        Rdiag = diag(R);
        Rdiag(end+1:end+3) = [(fOpt.sigmaUwbMPLA)^2 ...
                              (fOpt.sigmaUwbMPRA)^2 ...
                              (fOpt.sigmaUwbLARA)^2];
        R = diag(Rdiag);
    end
    
    alphalimit = struct('lkmin', fOpt.alphaLKmin, 'lkmax', fOpt.alphaLKmax, ...
                        'rkmin', fOpt.alphaRKmin, 'rkmax', fOpt.alphaRKmax);
    if (fOpt.applyCstr >= 1 && fOpt.applyCstr <= 8) || ...
        (fOpt.applyCstr >= 71 && fOpt.applyCstr <= 78) || ...
        (fOpt.applyCstr >= 201 && fOpt.applyCstr <= 208) || ...
        (fOpt.applyCstr >= 271 && fOpt.applyCstr <= 278)
        floorZ = min([x0(idxPosLA(3)), x0(idxPosRA(3))]);
    
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
%         if fOpt.applyAccBias
%             F(idxPosMP,idxAccBiasMP) = -dt2.*PELV_CS;
%             F(idxPosLA,idxAccBiasLA) = -dt2.*LTIB_CS;
%             F(idxPosRA,idxAccBiasRA) = -dt2.*RTIB_CS;
%             F(idxVelMP,idxAccBiasMP) = -dt.*PELV_CS;
%             F(idxVelLA,idxAccBiasLA) = -dt.*LTIB_CS;
%             F(idxVelRA,idxAccBiasRA) = -dt.*RTIB_CS;
%         end
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
        if mod(fOpt.applyMeas, 100) >= 1
            if bIsStatMP(n) idx(end+1:end+3) = idxMVelMP; end
            if bIsStatLA(n) idx(end+1:end+3) = idxMVelLA; end
            if bIsStatRA(n) idx(end+1:end+3) = idxMVelRA; end
        end
        
        if (mod(fOpt.applyMeas, 100) == 2) || ...
           (mod(fOpt.applyMeas, 100) >= 31 && mod(fOpt.applyMeas, 100) <= 33) || ...
           (mod(fOpt.applyMeas, 100) >= 41 && mod(fOpt.applyMeas, 100) <= 43) || ...
           (mod(fOpt.applyMeas, 100) >= 51 && mod(fOpt.applyMeas, 100) <= 53)
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
            end
        elseif mod(fOpt.applyMeas, 100) == 3
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ 
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA;
                y_k(idxMPosLA(1:2), n) = x_min(idxPosLA(1:2));
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA;
                y_k(idxMPosRA(1:2), n) = x_min(idxPosRA(1:2));
            end
        elseif mod(fOpt.applyMeas, 100) == 4
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+1) = idxMPosLA(3);
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+1) = idxMPosRA(3); 
            end
        elseif mod(fOpt.applyMeas, 100) == 5
            refSide = 'N';
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
            end
            if refSide == 'N'
                if bIsStatLA(n), refSide = 'L';
                elseif bIsStatRA(n), refSide = 'R'; end
            elseif ~(bIsStatLA(n) | bIsStatRA(n))
                refSide = 'N';
            end
        elseif mod(fOpt.applyMeas, 100) == 6
            idx(end+1:end+length(idxMPosMP)) = idxMPosMP; 
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
            end
        elseif mod(fOpt.applyMeas, 100) == 7
            idx(end+1:end+length(idxMPosMP)) = idxMPosMP;
            
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA; 
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA; 
            end
            if refSide == 'N'
                if bIsStatLA(n), refSide = 'L';
                elseif bIsStatRA(n), refSide = 'R'; end
            elseif ~(bIsStatLA(n) | bIsStatRA(n))
                refSide = 'N';
            end
        elseif mod(fOpt.applyMeas, 100) == 8
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ 
                idx(end+1:end+length(idxMPosLA)) = idxMPosLA;
                y_k(idxMPosLA(1:2), n) = x_min(idxPosLA(1:2));
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+length(idxMPosRA)) = idxMPosRA;
                y_k(idxMPosRA(1:2), n) = x_min(idxPosRA(1:2));
            end
        elseif mod(fOpt.applyMeas, 100) == 21
            if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                idx(end+1:end+1) = idxMPosLA;
            end
            if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                idx(end+1:end+1) = idxMPosRA; 
            end
            idx(end+1:end+6) = [idxMLHjc idxMRHjc];
            
            LTIB_CS = quat2rotm(y_k(idxMOriLA, n)');
            RTIB_CS = quat2rotm(y_k(idxMOriRA, n)');
            PELV_CS = quat2rotm(y_k(idxMOriMP, n)');
        
            y_k([idxMLHjc idxMRHjc], n) = solve_linhjc_d(x_min(idxPosMP,1), ...
                x_min(idxPosLA,1), x_min(idxPosRA,1), ...
                PELV_CS, LTIB_CS, RTIB_CS, ...
                dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
        end
        
        if (mod(fOpt.applyMeas, 100) >= 51 && mod(fOpt.applyMeas, 100) <= 53)
            modTenApplyMeas = mod(fOpt.applyMeas, 10);
            if bitand(modTenApplyMeas, 1)
                idx(end+1:end+3) = idxMVelMPLARA;
            end
            if bitand(modTenApplyMeas, 2)
                idx(end+1:end+2) = idxMPosMPLARA;
            end
        end
        
        if fOpt.applyMeas >= 100 && fOpt.applyMeas < 200
            idx(end+1:end+3) = [idxM3Dist];
            
            diff_MP_LA = x_min(idxPosMP)'-x_min(idxPosLA)';
            diff_MP_RA = x_min(idxPosMP)'-x_min(idxPosRA)';
            diff_LA_RA = x_min(idxPosLA)'-x_min(idxPosRA)';
            % the observation model
            h_uwb_est = [vecnorm(diff_MP_LA, 2, 2);
                         vecnorm(diff_MP_RA, 2, 2);
                         vecnorm(diff_LA_RA, 2, 2)];

            Huwb = zeros(3, nStates);
            % Jacobian of observation model with respect to elements of the state
            % estimate xhat, the order of these matrix elements MUST be preserved
            Huwb(1,1) = diff_MP_LA(1)/h_uwb_est(1);
            Huwb(1,2) = diff_MP_LA(2)/h_uwb_est(1);
            Huwb(1,3) = diff_MP_LA(3)/h_uwb_est(1);
            Huwb(1,7) = -diff_MP_LA(1)/h_uwb_est(1);
            Huwb(1,8) = -diff_MP_LA(2)/h_uwb_est(1);
            Huwb(1,9) = -diff_MP_LA(3)/h_uwb_est(1);    

            Huwb(2,1)  =  diff_MP_RA(1)/h_uwb_est(2);
            Huwb(2,2)  =  diff_MP_RA(2)/h_uwb_est(2);
            Huwb(2,3)  =  diff_MP_RA(3)/h_uwb_est(2);
            Huwb(2,13) = -diff_MP_RA(1)/h_uwb_est(2);
            Huwb(2,14) = -diff_MP_RA(2)/h_uwb_est(2);
            Huwb(2,15) = -diff_MP_RA(3)/h_uwb_est(2);

            Huwb(3,7)  =  diff_LA_RA(1)/h_uwb_est(3);
            Huwb(3,8)  =  diff_LA_RA(2)/h_uwb_est(3);
            Huwb(3,9)  =  diff_LA_RA(3)/h_uwb_est(3);
            Huwb(3,13) = -diff_LA_RA(1)/h_uwb_est(3);
            Huwb(3,14) = -diff_LA_RA(2)/h_uwb_est(3);
            Huwb(3,15) = -diff_LA_RA(3)/h_uwb_est(3);
            
            H(idxM3Dist, :) = Huwb;
            
            % calculate the measurement residual, i.e., the difference between the
            % predicted measurements from the observation model, h(xhat), and the
            % measurements obtained directly from the UWB sensors.
            y_k(idxM3Dist, n) = uwb_MP_LA_RA(n,:)' - h_uwb_est + Huwb * x_min;
        end
       
        res = y_k(idx, n) - H(idx, :) * x_min;
        K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx));
        x_min1 = x_min + K * res;
        
        if mod(fOpt.applyMeas,100) == 31
            if bIsStatLA(n), idx2 = [idx idxMCov1];
            else, idx2 = idx; end
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        elseif mod(fOpt.applyMeas,100) == 32
            idx2 = [idx idxMCov1];
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        elseif mod(fOpt.applyMeas,100) == 33
            if bIsStatLA(n), idx2 = [idx idxMCov2];
            else, idx2 = [idx idxMCov1]; end
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
        else
            idx2 = idx;
        end
        P_min1 = (I_N - K * H(idx2, :)) * P_min;
        
        if mod(fOpt.applyMeas, 100) == 4
            idx = [idxMPosLA(1:2) idxMPosRA(1:2)];
            y_k(idx) = [x_min(idxPosLA(1:2)) x_min(idxPosRA(1:2))];
            K = P_min * H(idx, :)' /(H(idx, :) * P_min * H(idx,:)' + R(idx, idx));

            P_min1 = (I_N - K * H(idx, :)) * P_min1;
        elseif mod(fOpt.applyMeas, 100) == 5 || mod(fOpt.applyMeas, 100) == 7
            % zero out part of P_min1
            if refSide == 'L'
                P_min1(idxPosLA, :) = 0;
                P_min1(:, idxPosLA) = 0;
            elseif refSide == 'R'
                P_min1(idxPosRA, :) = 0;
                P_min1(:, idxPosRA) = 0;
            end
        elseif mod(fOpt.applyMeas,100) == 41
            if bIsStatLA(n)
                idx2 = idxMCov1;
                K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
                P_min1 = (I_N - K * H(idx2, :)) * P_min1;
            end
        elseif mod(fOpt.applyMeas,100) == 42
            idx2 = idxMCov1;
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
            P_min1 = (I_N - K * H(idx2, :)) * P_min1;
        elseif mod(fOpt.applyMeas,100) == 43
            if bIsStatLA(n), idx2 = idxMCov2;
            else, idx2 = idxMCov1; end
            K = P_min * H(idx2, :)' /(H(idx2, :) * P_min * H(idx2,:)' + R(idx2, idx2));
            P_min1 = (I_N - K * H(idx2, :)) * P_min1;
        end
        
        if fOpt.applyMeas
            debug_dat.zuptState(n,:) = x_min1;
            debug_dat.zuptP(:,:,n) = P_min1;
        end
        
    %% ---- Kalman Filter Update Step using UWB measurements ---- 
    % this correction step should be done last
        if fOpt.applyUwb
            diff_MP_LA = x_min1(idxPosMP)' - x_min1(idxPosLA)';
            diff_MP_RA = x_min1(idxPosMP)' - x_min1(idxPosRA)';
            diff_LA_RA = x_min1(idxPosLA)' - x_min1(idxPosRA)';
            % the observation model
            h_uwb_est = [vecnormalize( diff_MP_LA );
                         vecnormalize( diff_MP_RA );
                         vecnormalize( diff_LA_RA );];
            % calculate the measurement residual, i.e., the difference between the
            % predicted measurements from the observation model, h(xhat), and the
            % measurements obtained directly from the UWB sensors.
            y_innovation = uwb_MP_LA_RA(n,:)' - h_uwb_est;

            H_uwb = zeros(3, nStates);
            % Jacobian of observation model with respect to elements of the state
            % estimate xhat, the order of these matrix elements MUST be preserved
            H_uwb(1,1) = diff_MP_LA(1)/h_uwb_est(1);
            H_uwb(1,2) = diff_MP_LA(2)/h_uwb_est(1);
            H_uwb(1,3) = diff_MP_LA(3)/h_uwb_est(1);
            H_uwb(1,7) = -diff_MP_LA(1)/h_uwb_est(1);
            H_uwb(1,8) = -diff_MP_LA(2)/h_uwb_est(1);
            H_uwb(1,9) = -diff_MP_LA(3)/h_uwb_est(1);    

            H_uwb(2,1)  =  diff_MP_RA(1)/h_uwb_est(2);
            H_uwb(2,2)  =  diff_MP_RA(2)/h_uwb_est(2);
            H_uwb(2,3)  =  diff_MP_RA(3)/h_uwb_est(2);
            H_uwb(2,13) = -diff_MP_RA(1)/h_uwb_est(2);
            H_uwb(2,14) = -diff_MP_RA(2)/h_uwb_est(2);
            H_uwb(2,15) = -diff_MP_RA(3)/h_uwb_est(2);

            H_uwb(3,7)  =  diff_LA_RA(1)/h_uwb_est(3);
            H_uwb(3,8)  =  diff_LA_RA(2)/h_uwb_est(3);
            H_uwb(3,9)  =  diff_LA_RA(3)/h_uwb_est(3);
            H_uwb(3,13) = -diff_LA_RA(1)/h_uwb_est(3);
            H_uwb(3,14) = -diff_LA_RA(2)/h_uwb_est(3);
            H_uwb(3,15) = -diff_LA_RA(3)/h_uwb_est(3);    
            % Calculate the covariance in the measurement residual
            S    = ((H_uwb * P_min1) * H_uwb') + R_uwb; % scalar
            % Calculate the kalman gain
            K    = P_min1 * H_uwb' * S^(-1);
            % Update state estimate
            x_plus = x_min1 + K * y_innovation;
            % Update Covariance in the state estimate
            P_plus    = (I_N - K * H_uwb) * P_min1;

            debug_dat.uwbuptState(n,:) = x_plus;
            debug_dat.uwbuptP(:,:,n) = P_plus;
        else
            x_plus = x_min1;
            P_plus = P_min1;
        end
        
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
        
        if fOpt.applyCstr == 0
            x_tilde = x_plus;
            P_tilde = P_plus;
        elseif (fOpt.applyCstr >= 1 && fOpt.applyCstr <= 8 ) || ...
            (fOpt.applyCstr >= 71 && fOpt.applyCstr <= 78 ) || ...
            (fOpt.applyCstr >= 201 && fOpt.applyCstr <= 208 ) || ...
            (fOpt.applyCstr >= 271 && fOpt.applyCstr <= 278 )
            if fOpt.applyCstr < 100
                d_k = solve_linhjc_d(x_plus(idxPosMP,1), x_plus(idxPosLA,1), ...
                        x_plus(idxPosRA,1), PELV_CS, LTIB_CS, RTIB_CS, ...
                        dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
            else
                d_k = solve_linhjc_kac_d(x_plus(idxPosMP,1), x_plus(idxPosLA,1), ...
                        x_plus(idxPosRA,1), PELV_CS, LTIB_CS, RTIB_CS, ...
                        dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
            end
            
            switch (mod(fOpt.applyCstr, 10))
                case 1 % maximum probability estimate + no P update
                    Kk = P_plus*D'*(D*P_plus*D')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = P_plus;
                case 2 % least squares estimate + no P update
                    Kk = D'*(D*D')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = P_plus;
                case 3 % least squares estimate w/ full confidence on pelvis (P=0 at pelvis) + no P update
                    Kk = P_custom*D'*(D*P_custom*D')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = P_plus;
                case 4 % maximum probability estimate w/ force equal foot covariance + no P update
                    Winv = P_plus;
                    Winv(idxPosLA, idxPosLA) = Winv(idxPosRA, idxPosRA);
                    Winv(idxVelLA, idxVelLA) = Winv(idxVelRA, idxVelRA);
                    Kk = Winv*D'*(D*Winv*D')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = P_plus;
                case 5 % maximum probability estimate if P is well conditioned + P update
                    if cond(P_plus) < 1e5
                        Kk = P_plus*D'*(D*P_plus*D')^(-1);
                        x_tilde = x_plus + Kk*(d_k - D * x_plus);
                        P_tilde = (I_N-Kk*D)*P_plus*(I_N-Kk*D)';
                    else
                        x_tilde = x_plus;
                        P_tilde = P_plus;
                        debug_dat.cstrStateU(n, 1) = false;
                    end
                case 6 % soft maximum probability estimate + P update
                    Kk = P_plus * D' * (D * P_plus * D' + fOpt.sigmaCPos*I_CN)^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = (I_N-Kk*D)*P_plus*(I_N-Kk*D)' + Kk*fOpt.sigmaCPos*I_CN*Kk';
                case 7 % maximum probability estimate of constraint subset + P update
                    P_cstr = D*P_plus*D';
                    idx = diag(P_cstr) > 1e-5;
                    
                    if sum(idx) > 0
                        D2 = D(idx, :);
                        Kk = P_plus*D2'*(P_cstr(idx, idx))^(-1);
                        x_tilde = x_plus + Kk*(d_k(idx,:) - D2 * x_plus);
                        P_tilde = (I_N-Kk*D2)*P_plus*(I_N-Kk*D2)';
                    else
                        x_tilde = x_plus;
                        P_tilde = P_plus;
                    end
                    debug_dat.cstrStateU(n, :) = idx;
                case 8 % maximum probability estimate + soft P update
                    Kk = P_plus * D' * (D * P_plus * D')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D * x_plus);
                    P_tilde = (I_N-Kk*D)*P_plus*(I_N-Kk*D)' + Kk*fOpt.sigmaCPos*I_CN*Kk';
            end
        elseif fOpt.applyCstr >= 11 && fOpt.applyCstr <= 16 % fmincon linear const
            % calculate the location of the knee
            
            d_k = solve_linhjc_d(x_plus(idxPosMP,1), x_plus(idxPosLA,1), ...
                    x_plus(idxPosRA,1), PELV_CS, LTIB_CS, RTIB_CS, ...
                    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
            x_hat = x_plus(x2cxIdx,1);
            x_hat0 = x_plus(x2cxIdx,1);
            
            switch (fOpt.applyCstr)
                case 11 % fmincon (interior point) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'interior-point';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 12 % fmincon (sqp) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'sqp';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 13 % fmincon (active-set) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'active-set';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 14 % fmincon (interior point) linear hjc, W = I
                    optimOpt.Algorithm = 'interior-point';
                    P_hat = eye(nCStates);
                case 15 % fmincon (sqp) linear hjc, W = I
                    optimOpt.Algorithm = 'sqp';
                    P_hat = eye(nCStates);
                case 16 % fmincon (active-set) linear hjc, W = I
                    optimOpt.Algorithm = 'active-set';
                    P_hat = eye(nCStates);
            end
            
            x_hat = fmincon(@(x_hat) L2(x_hat, x_hat0, P_hat, 20), ...
                            x_hat0, [], [], D, d_k, [], [], [], optimOpt);
                            
            x_tilde = x_plus;
            x_tilde(x2cxIdx, 1) = x_hat;            
            P_tilde = P_plus;
        elseif (fOpt.applyCstr >= 21 && fOpt.applyCstr <= 23) || ...
               (fOpt.applyCstr >= 221 && fOpt.applyCstr <= 223)
            x_plus2 = pelib.est.changeStateRefFrame(x_plus', 'MIDPEL', ...
                            'kf_3_kmus_v3')';
            
            LTIB_CS2 = quat2rotm(x_plus2(idxOriLA,1)');
            RTIB_CS2 = quat2rotm(x_plus2(idxOriRA,1)');
            PELV_CS2 = quat2rotm(x_plus2(idxOriMP,1)');
        
            if fOpt.applyCstr < 100
                d_k = solve_linhjc_d(x_plus2(idxPosMP,1), x_plus2(idxPosLA,1), ...
                    x_plus2(idxPosRA,1), PELV_CS2, LTIB_CS2, RTIB_CS2, ...
                    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
            else
                d_k = solve_linhjc_kac_d(x_plus2(idxPosMP,1), x_plus2(idxPosLA,1), ...
                    x_plus2(idxPosRA,1), PELV_CS2, LTIB_CS2, RTIB_CS2, ...
                    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
            end           
            
            switch (mod(fOpt.applyCstr, 100))
                case 21 % maximum probability estimate + no P update
                    Kk = P_plus*D'*(D*P_plus*D')^(-1);
                    x_tilde2 = x_plus2 + Kk*(d_k - D * x_plus2);
                    P_tilde = P_plus;
                case 22 % least squares estimate + no P update
                    Kk = D'*(D*D')^(-1);
                    x_tilde2 = x_plus2 + Kk*(d_k - D * x_plus2);
                    P_tilde = P_plus;
                case 23 % maximum probability estimate of constraint subset + P update
                    P_cstr = D*P_plus*D';
                    idx = diag(P_cstr) > 1e-5;
                    
                    if sum(idx) > 0
                        D2 = D(idx, :);
                        Kk = P_plus*D2'*(P_cstr(idx, idx))^(-1);
                        x_tilde2 = x_plus2 + Kk*(d_k(idx,:) - D2 * x_plus2);
                        P_tilde = (I_N-Kk*D2)*P_plus*(I_N-Kk*D2)';
                    else
                        x_tilde2 = x_plus2;
                        P_tilde = P_plus;
                    end
                    debug_dat.cstrStateU(n, :) = idx;
            end
            
            x_tilde = x_plus;
            x_tilde(idxPosLA,1) = x_plus(idxPosMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_tilde2(idxPosLA,1)')';
            x_tilde(idxVelLA,1) = x_plus(idxVelMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_tilde2(idxVelLA,1)')';
            x_tilde(idxPosRA,1) = x_plus(idxPosMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_tilde2(idxPosRA,1)')';
            x_tilde(idxVelRA,1) = x_plus(idxVelMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_tilde2(idxVelRA,1)')';
        elseif fOpt.applyCstr >= 31 && fOpt.applyCstr <= 36 % fmincon linear const
            % calculate the location of the knee
            
            x_plus2 = pelib.est.changeStateRefFrame(x_plus', 'MIDPEL', ...
                            'kf_3_kmus_v3')';
            
            LTIB_CS2 = quat2rotm(x_plus2(idxOriLA,1)');
            RTIB_CS2 = quat2rotm(x_plus2(idxOriRA,1)');
            PELV_CS2 = quat2rotm(x_plus2(idxOriMP,1)');
        
            d_k = solve_linhjc_d(x_plus2(idxPosMP,1), x_plus2(idxPosLA,1), ...
                    x_plus2(idxPosRA,1), PELV_CS2, LTIB_CS2, RTIB_CS2, ...
                    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
                
            x_hat0 = x_plus2(x2cxIdx,1);
            
            switch (fOpt.applyCstr)
                case 31 % fmincon (interior point) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'interior-point';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 32 % fmincon (sqp) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'sqp';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 33 % fmincon (active-set) linear hjc, W = P^-1
                    optimOpt.Algorithm = 'active-set';
                    P_hat = P_plus(x2cxIdx, x2cxIdx);
                case 34 % fmincon (interior point) linear hjc, W = I
                    optimOpt.Algorithm = 'interior-point';
                    P_hat = eye(nCStates);
                case 35 % fmincon (sqp) linear hjc, W = I
                    optimOpt.Algorithm = 'sqp';
                    P_hat = eye(nCStates);
                case 36 % fmincon (active-set) linear hjc, W = I
                    optimOpt.Algorithm = 'active-set';
                    P_hat = eye(nCStates);
            end
            
            x_hat = fmincon(@(x_hat) L2(x_hat, x_hat0, P_hat, 20), ...
                            x_hat0, [], [], D, d_k, [], [], [], optimOpt);
                            
            x_tilde = x_plus;
            x_tilde(idxPosLA,1) = x_plus(idxPosMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_hat(1:3,1)')';
            x_tilde(idxVelLA,1) = x_plus(idxVelMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_hat(4:6,1)')';
            x_tilde(idxPosRA,1) = x_plus(idxPosMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_hat(7:9,1)')';
            x_tilde(idxVelRA,1) = x_plus(idxVelMP,1) + ...
                quatrotate(quatconj(x_plus(idxOriMP,1)'), x_hat(10:12,1)')';
            
            P_tilde = P_plus;
        elseif fOpt.applyCstr >= 141 && fOpt.applyCstr <= 146
            sckfAlpha = fOpt.sckfAlpha;
            sckfThreshold = fOpt.sckfThreshold;
            x_tilde = x_plus;
            P_tilde = P_plus;
            
            for i=0:fOpt.sckfMaxIter
                % calculate the z axis of femur and tibia
                LFEM_z = x_tilde(idxPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxPosLA,1);
                RFEM_z = x_tilde(idxPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxPosRA,1);

                g_dlfem = norm(LFEM_z, 2);
                g_drfem = norm(RFEM_z, 2);

                D = zeros(4, length(x_tilde));
                D(1, idxPosMP) = LFEM_z'/g_dlfem;
                D(1, idxPosLA) = -LFEM_z'/g_dlfem;
                D(2, idxPosMP) = RFEM_z'/g_drfem;
                D(2, idxPosRA) = -RFEM_z'/g_drfem;
                D(3, idxPosMP) = LTIB_CS(:,2)';
                D(3, idxPosLA) = -LTIB_CS(:,2)';
                D(4, idxPosMP) = RTIB_CS(:,2)';
                D(4, idxPosRA) = -RTIB_CS(:,2)';

                if i==0
                    R0 = sckfAlpha*D*P_tilde*D';
                end

                Ri = R0*exp(-i);
                Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ diag(D * P_tilde* D');
                Si(isnan(Si)) = sckfThreshold+1;
                if sum(Si < sckfThreshold) == 0, break, end

                switch (fOpt.applyCstr)
                    case 141
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 142
                        Kk = D'*(D*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 143
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                    case 144
                        Kk = D'*(D*D'+Ri)^(-1);
                    case 145
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 146
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                    otherwise
                        Kk = 0;
                end
                
                res = [dLFemur - g_dlfem; dRFemur - g_drfem; ...
                       (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*LTIB_CS(:,2) ...
                            - D(3,:)*x_tilde; ...
                       (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*RTIB_CS(:,2) ...
                            - D(4,:)*x_tilde];
                dx = Kk*(res);
                x_tilde = x_tilde + dx;
            end
            P_tilde = P_plus;
        elseif (fOpt.applyCstr >= 151 && fOpt.applyCstr <= 158) || ...
               (fOpt.applyCstr >= 161 && fOpt.applyCstr <= 168) || ...
               (fOpt.applyCstr >= 171 && fOpt.applyCstr <= 178)
            sckfAlpha = fOpt.sckfAlpha;
            sckfThreshold = fOpt.sckfThreshold;
            x_tilde = x_plus;
            P_tilde = P_plus;
            
            for i=0:fOpt.sckfMaxIter
                % calculate the z axis of femur and tibia
                LFEM_z = x_tilde(idxPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxPosLA,1);
                RFEM_z = x_tilde(idxPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxPosRA,1);
                if i==0
                    alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                                 -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
                    alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                                     -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
                                 
                    if mod(fOpt.applyCstr, 10) >= 7 && mod(fOpt.applyCstr, 10) <= 8 
%                         alphalimit.lkmin = max(alpha_lk, fOpt.alphaLKmin);
%                         alphalimit.rkmin = max(alpha_rk, fOpt.alphaRKmin);
                        
                        if bIsStatLA(n) || x_min(idxPosLA(3), 1) < floorZ
                            P_tilde(idxPosLA(3), :) = 0; 
                            P_tilde(:, idxPosLA(3)) = 0; 
                        end
                        if bIsStatRA(n) || x_min(idxPosRA(3), 1) < floorZ
                            P_tilde(idxPosRA(3), :) = 0; 
                            P_tilde(:, idxPosRA(3)) = 0; 
                        end
                    end
                    
                    if alpha_lk < alphalimit.lkmin
                        tmpLK = sin(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,1) - ...
                                cos(alphalimit.lkmin - 0.5*pi)*LTIB_CS(:,3);
                    elseif alpha_lk > alphalimit.lkmax
                        tmpLK = sin(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,1) - ...
                                cos(alphalimit.lkmax - 0.5*pi)*LTIB_CS(:,3);
                    else
                        tmpLK = false;
                    end
                    
                    if alpha_rk < alphalimit.rkmin
                        tmpRK = sin(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,1) - ...
                                cos(alphalimit.rkmin - 0.5*pi)*RTIB_CS(:,3);
                    elseif alpha_rk > alphalimit.rkmax
                        tmpRK = sin(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,1) - ...
                                cos(alphalimit.rkmax - 0.5*pi)*RTIB_CS(:,3);
                    else
                        tmpRK = false;
                    end
                end
                
                g_dlfem = norm(LFEM_z, 2);
                g_drfem = norm(RFEM_z, 2);

                D = zeros(4, length(x_tilde));
                D(1, idxPosMP) = LFEM_z'/g_dlfem;
                D(1, idxPosLA) = -LFEM_z'/g_dlfem;
                D(2, idxPosMP) = RFEM_z'/g_drfem;
                D(2, idxPosRA) = -RFEM_z'/g_drfem;
                D(3, idxPosMP) = LTIB_CS(:,2)';
                D(3, idxPosLA) = -LTIB_CS(:,2)';
                D(4, idxPosMP) = RTIB_CS(:,2)';
                D(4, idxPosRA) = -RTIB_CS(:,2)';
                
                res = [dLFemur - g_dlfem; dRFemur - g_drfem; ...
                       (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*LTIB_CS(:,2) ...
                            - D(3,:)*x_tilde; ...
                       (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*RTIB_CS(:,2) ...
                            - D(4,:)*x_tilde];
                        
                if isnumeric(tmpLK)
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpLK';
                    D(end, idxPosLA) = -tmpLK';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLK ...
                            - D(end,:)*x_tilde;
                end
                
                if isnumeric(tmpRK)
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpRK';
                    D(end, idxPosRA) = -tmpRK';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRK ...
                            - D(end,:)*x_tilde;
                end
                
                if i==0
                    R0 = sckfAlpha*D*P_tilde*D';
                end

                Ri = R0*exp(-i);
                Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ diag(D * P_tilde* D');
                Si(isnan(Si)) = sckfThreshold+1;
                if sum(Si < sckfThreshold) == 0, break, end

                switch mod(fOpt.applyCstr, 10)
                    case 1
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 2
                        Kk = D'*(D*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 3
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                    case 4
                        Kk = D'*(D*D'+Ri)^(-1);
                    case 5
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 6
                        Kk = P_custom*D'*(D*P_custom*D'+Ri)^(-1);
                    case 7
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 8
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                    otherwise
                        Kk = 0;
                end
                
                dx = Kk*(res);
                x_tilde = x_tilde + dx;
            end
            if ~(fOpt.applyCstr >= 161 && fOpt.applyCstr <= 168)
                P_tilde = P_plus;
            end
        elseif (fOpt.applyCstr >= 351 && fOpt.applyCstr <= 356) || ...
               (fOpt.applyCstr >= 371 && fOpt.applyCstr <= 376) || ...
               (fOpt.applyCstr >= 381 && fOpt.applyCstr <= 386)
            sckfAlpha = fOpt.sckfAlpha;
            sckfThreshold = fOpt.sckfThreshold;
            x_tilde = x_plus;
            P_tilde = P_plus;
            
            % preprocessing for knee angle inequality
            % additional process to ensure knee angle does not increase
            % during iteration
            if (mod(fOpt.applyCstr, 10) >= 5 && mod(fOpt.applyCstr, 10) <= 6)
                LFEM_z = x_tilde(idxPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxPosLA,1);
                RFEM_z = x_tilde(idxPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxPosRA,1);
                     
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
                LFEM_z = x_tilde(idxPosMP,1) + dPelvis/2*PELV_CS(:,2) ...
                         - dLTibia*LTIB_CS(:,3) - x_tilde(idxPosLA,1);
                RFEM_z = x_tilde(idxPosMP,1) - dPelvis/2*PELV_CS(:,2) ...
                         - dRTibia*RTIB_CS(:,3) - x_tilde(idxPosRA,1);
                g_dlfem = norm(LFEM_z, 2);
                g_drfem = norm(RFEM_z, 2);
                
                % D matrix construction: femur length constraint and knee hinge joint
                D = zeros(4, length(x_tilde));
                D(1, idxPosMP) = LFEM_z'/g_dlfem;
                D(1, idxPosLA) = -LFEM_z'/g_dlfem;
                D(2, idxPosMP) = RFEM_z'/g_drfem;
                D(2, idxPosRA) = -RFEM_z'/g_drfem;
                D(3, idxPosMP) = LTIB_CS(:,2)';
                D(3, idxPosLA) = -LTIB_CS(:,2)';
                D(4, idxPosMP) = RTIB_CS(:,2)';
                D(4, idxPosRA) = -RTIB_CS(:,2)';

                res = [dLFemur - g_dlfem; ...
                       dRFemur - g_drfem; ...
                      (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*LTIB_CS(:,2) ...
                            - D(3,:)*x_tilde; ...
                      (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*RTIB_CS(:,2) ...
                            - D(4,:)*x_tilde];
                  
                % add knee inequality constraint via active set
                alpha_lk = atan2(-dot(LFEM_z, LTIB_CS(:,3)), ...
                                 -dot(LFEM_z, LTIB_CS(:,1))) + 0.5*pi;
                alpha_rk = atan2(-dot(RFEM_z, RTIB_CS(:,3)), ...
                                 -dot(RFEM_z, RTIB_CS(:,1))) + 0.5*pi;
                if alpha_lk < alphalimit.lkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpLKmin';
                    D(end, idxPosLA) = -tmpLKmin';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmin ...
                        - D(end, :)*x_tilde;
                elseif alpha_lk > alphalimit.lkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpLKmax';
                    D(end, idxPosLA) = -tmpLKmax';
                    res(end+1, 1) = (-dPelvis/2*PELV_CS(:,2) + dLTibia*LTIB_CS(:,3))'*tmpLKmax ...
                        - D(end, :)*x_tilde;
                end
                    
                if alpha_rk < alphalimit.rkmin
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpRKmin';
                    D(end, idxPosRA) = -tmpRKmin';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmin ...
                        - D(end, :)*x_tilde;
                elseif alpha_rk > alphalimit.rkmax
                    D(end+1, :) = zeros(1, length(x_tilde));
                    D(end, idxPosMP) = tmpRKmax';
                    D(end, idxPosRA) = -tmpRKmax';
                    res(end+1, 1) = (dPelvis/2*PELV_CS(:,2) + dRTibia*RTIB_CS(:,3))'*tmpRKmax ...
                        - D(end, :)*x_tilde;
                end
                
                % add flat floor constraint
                if (fOpt.applyCstr >= 371 && fOpt.applyCstr <= 376)
                    stepCstrCount = 0;
                    if bIsStatLA(n) || x_tilde(idxPosLA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxPosLA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        stepCstrCount = stepCstrCount + 1;
                    end
                    if bIsStatRA(n) || x_tilde(idxPosRA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxPosRA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        stepCstrCount = stepCstrCount + 1;
                    end
                    if stepCstrCount == 2
                        P_tilde_cstr = D*P_tilde*D';
                        % if condition number is high enough (singular
                        % matrix), delete last constraint
                        if cond(P_tilde_cstr) > 1e5
                            D(end, :) = [];
                            res(end, :) = [];
                        end
                    end
                elseif(fOpt.applyCstr >= 381 && fOpt.applyCstr <= 386)
                    stepCstrCount = 0;
                    if bIsStatLA(n) || x_tilde(idxPosLA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxPosLA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        stepCstrCount = stepCstrCount + 1;
                    end
                    if bIsStatRA(n) || x_tilde(idxPosRA(3), 1) < floorZ
                        D(end+1, :) = zeros(1, length(x_tilde));
                        D(end, idxPosRA(3)) = 1;
                        res(end+1, 1) = floorZ - D(end, :)*x_tilde;
                        stepCstrCount = stepCstrCount + 1;
                    end
                    if stepCstrCount == 2
                        P_tilde_cstr = D*P_tilde*D';
                        % if condition number is high enough (singular
                        % matrix), delete last constraint
                        if cond(P_tilde_cstr) > 1e5
                            D(end, :) = [];
                            res(end, :) = [];
                        end
                    end
                    if bIsStatLA(n)
                        D(end+1:end+2, :) = zeros(2, length(x_tilde));
                        D(end-1:end, idxPosLA(1:2)) = eye(2,2);
                        res(end+1:end+2, 1) = x_plus(idxPosLA(1:2), 1) - x_tilde(idxPosLA(1:2), 1);
                    elseif bIsStatRA(n)
                        D(end+1:end+2, :) = zeros(2, length(x_tilde));
                        D(end-1:end, idxPosRA(1:2)) = eye(2,2);
                        res(end+1:end+2, 1) = x_plus(idxPosRA(1:2), 1) - x_tilde(idxPosRA(1:2), 1);
                    end
                end
                
                % Step 3: calculate Ri               
                Ri = sckfAlpha*D*P_tilde*D'*exp(-i);
                Si_denom = diag(D * P_tilde* D');
                Si = max(D.^2 .* diag(P_tilde)', [], 2) ./ Si_denom;
                Si_ignoreindex = isnan(Si) | (Si_denom < 1e-10);
                Si(Si_ignoreindex) = sckfThreshold+1;
                if sum(Si < sckfThreshold) == 0, break, end

                % Step 4: update x_tilde and P_tilde
                switch mod(fOpt.applyCstr, 10)
                    case 1
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 2
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                    case 3
                        Kk = D'*(D*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 4
                        Kk = D'*(D*D'+Ri)^(-1);
                    case 5
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                        P_tilde = (I_N-Kk*D)*P_tilde*(I_N-Kk*D)' + Kk*Ri*Kk';
                    case 6
                        Kk = P_tilde*D'*(D*P_tilde*D'+Ri)^(-1);
                    otherwise
                        Kk = 0;
                end
                
                dx = Kk*(res);
                x_tilde = x_tilde + dx;
            end
            P_tilde = P_plus;
        elseif (fOpt.applyCstr >= 51 && fOpt.applyCstr <= 54)
            % 001 constraints + MP/LA/RA zpos = floor zpos         
            d_k = solve_linhjc_d(x_plus(idxPosMP,1), x_plus(idxPosLA,1), ...
                    x_plus(idxPosRA,1), PELV_CS, LTIB_CS, RTIB_CS, ...
                    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia);
                
            idx = 1:6;
            lowestpoints = x_plus([idxPosMP(3), idxPosLA(3), idxPosRA(3)]);
            [lpy lpi] = min(lowestpoints);
            if lpy < floorZ
                idx(end+1) = lpi+6;
                d_k(end+1) = floorZ;
            elseif lpi == 2 && bIsStatLA(n)
                idx(end+1) = 8;
                d_k(end+1) = floorZ;
            elseif lpi == 3 && bIsStatRA(n)
                idx(end+1) = 9;
                d_k(end+1) = floorZ;
            end
            
            D2 = D(idx,:);
            switch (fOpt.applyCstr)
                case 51 % maximum probability estimate + no P update
                    Kk = P_plus*D2'*(D2*P_plus*D2')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D2 * x_plus);
                    P_tilde = P_plus;
                case 52 % least squares estimate + no P update
                    Kk = D2'*(D2*D2')^(-1);
                    x_tilde = x_plus + Kk*(d_k - D2 * x_plus);
                    P_tilde = P_plus;
                case 53 % soft maximum probability estimate + P update
                    R_tilde = fOpt.sigmaCPos*eye(length(idx));
                    Kk = P_plus * D2' * (D2 * P_plus * D2' + R_tilde)^(-1);
                    x_tilde = x_plus + Kk*(d_k - D2 * x_plus);
                    P_tilde = (I_N-Kk*D2)*P_plus*(I_N-Kk*D2)' + Kk*R_tilde*Kk';
                case 54 % maximum probability estimate of constraint subset + P update
                    P_cstr = D2*P_plus*D2';
                    idx2 = diag(P_cstr) > 1e-5;
                    
                    if sum(idx2) > 0
                        D3 = D2(idx2, :);
                        Kk = P_plus*D3'*(P_cstr(idx2, idx2))^(-1);
                        x_tilde = x_plus + Kk*(d_k(idx2,:) - D3 * x_plus);
                        P_tilde = (I_N-Kk*D3)*P_plus*(I_N-Kk*D3)';
                    else
                        x_tilde = x_plus;
                        P_tilde = P_plus;
                    end
                    debug_dat.cstrStateU(n, 1:6) = idx2(1:6);
            end
            
            debug_dat.cstrStateU(n, idx(end)) = true;
        elseif (fOpt.applyCstr >= 101 && fOpt.applyCstr <= 106) || ... % fmincon
            (fOpt.applyCstr >= 111 && fOpt.applyCstr <= 116)
        
            switch (fOpt.applyCstr)
                case 101
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus;
                case 102
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus;
                case 103
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus;
                case 104
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = I_N;
                case 105
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = I_N;
                case 106
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x_tilde) hjc_nonlcon(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = I_N;
                case 111
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus;
                case 112
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus;
                case 113
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus;
                case 114
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = I_N;
                case 115
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = I_N;
                case 116
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x_tilde) hjc_nonlcon_kac(x_tilde, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = I_N;
            end
            
            x_tilde = fmincon(@(x_tilde) L2Dist(x_tilde, x_plus, W, 20), ...
                x_plus, [], [], [], [], [], [], nonlcon, optimOpt);
            
            % normalize the quaternions
            x_tilde(idxOriMP) = quatnormalize(x_tilde(idxOriMP, 1)');
            x_tilde(idxOriLA) = quatnormalize(x_tilde(idxOriLA, 1)');
            x_tilde(idxOriRA) = quatnormalize(x_tilde(idxOriRA, 1)');
            P_tilde = P_plus;
            
        elseif (fOpt.applyCstr >= 121 && fOpt.applyCstr <= 126) || ... % fmincon
            (fOpt.applyCstr >= 131 && fOpt.applyCstr <= 136)
        
            idx = [idxPosMP, idxPosLA, idxPosRA];
            switch (fOpt.applyCstr)
                case 121
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus(idx, idx);
                case 122
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus(idx, idx);
                case 123
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = P_plus(idx, idx);
                case 124
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = eye(9);
                case 125
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = eye(9);
                case 126
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x) hjc_nonlcon_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia);
                    W = eye(9);
                case 131
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus(idx, idx);
                case 132
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus(idx, idx);
                case 133
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = P_plus(idx, idx);
                case 134
                    optimOpt.Algorithm = 'interior-point';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = eye(9);
                case 135
                    optimOpt.Algorithm = 'sqp';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = eye(9);
                case 136
                    optimOpt.Algorithm = 'active-set';
                    nonlcon = @(x) hjc_nonlcon_kac_poska(x, ...
                                PELV_CS, LTIB_CS, RTIB_CS, dPelvis, ...
                                dLFemur, dRFemur, dLTibia, dRTibia, alphalimit);
                    W = eye(9);
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

        if (fOpt.applyCstr >= 71 && fOpt.applyCstr <= 77) || ...
           (fOpt.applyCstr >= 271 && fOpt.applyCstr <= 277) || ...
           (fOpt.applyCstr >= 171 && fOpt.applyCstr <= 178)
            % 001 constraints + MP/LA/RA zpos = floor zpos         
            idx = [idxPosMP(3), idxPosLA(3), idxPosRA(3)];
            [lpy lpi] = min(x_tilde(idx));
            if lpy < floorZ || (lpi == 2 && bIsStatLA(n)) || ...
                    (lpi == 3 && bIsStatRA(n))
                x_tilde(idx) = x_tilde(idx) - x_tilde(idx(lpi)) + floorZ;
            end
        end
        
        xhat_con(n, :) = x_tilde;
        P_con(:, :, n)  = P_tilde;
        debug_dat.cstrState(n,:) = x_tilde;
        debug_dat.cstrP(:,:,n) = P_tilde;
        
        debug_dat.LFEO(n, :) = x_tilde(idxPosLA) + dLTibia * LTIB_CS(:, 3);
        debug_dat.RFEO(n, :) = x_tilde(idxPosRA) + dRTibia * RTIB_CS(:, 3);
        debug_dat.LFEP(n, :) = x_tilde(idxPosMP) + dPelvis/2 * PELV_CS(:, 2);
        debug_dat.RFEP(n, :) = x_tilde(idxPosMP) - dPelvis/2 * PELV_CS(:, 2);
        
        if fOpt.applyCstr
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
        end
    end
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
                        dPelvis, dLFemur, dRFemur, dLTibia, dRTibia)
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