%> EKF using Lie group/algebra representation
%> 3 KMUs presumably worn on the body in the following configuration: 
%> mid pelvis, left ankle, right ankle
%>
%> Author: Luke Wicent Sy
%>
%> Inputs:
%>   x0       - the initial state in the GFR
%>   P0       - initial covariance
%>   gfrAccMP - the acceleration of the mid-pelvis in the GFR
%>   gfrAccLA - the acceleration of the left ankle in the GFR
%>   gfrAccRA - the acceleration of the right ankle in the GFR
%>   bIsStatMP  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_MP(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   bIsStatLA  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_LA(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   bIsStatRA  - a boolean vector, for whichever timepoints, n(i) are true,
%>                i.e., bMoving_RA(i) == 1, a zero velocity update will be 
%>                performed by using psuedo-zero velocity measurements 
%>   qMP       - mid  pelvis orientation in the GFR (quaternion)
%>   qLA       - left  ankle orientation in the GFR (quaternion)
%>   qRA       - right ankle orientation in the GFR (quaternion)
%>   wMP       - pelvis      angular velocity in the GFR
%>   wLA       - left  shank angular velocity in the GFR
%>   wRA       - right shank angular velocity in the GFR
%>   dPelvis   - pelvis width
%>   dRFemur   - right femur length
%>   dLFemur   - left femur length
%>   dRTibia   - right tibia length
%>   dLTibia   - left tibia length
%>   uwb_mea   - a structure containing the range measurements (m) between
%>   vel0      - struct of MP, LA, RA containing initial velocity of joints
%>   options   - struct containing the ff. settings:
%>   + fs      - sampling frequency of the magnetic and inertial measurement units

function [ xhat_pri, xhat_pos, debug_dat ] = lieekf_3_kmus_v3(x0, P0, ...
    gfrAccMP, bIsStatMP, qMP, wMP, ...
    gfrAccLA, bIsStatLA, qLA, wLA, ...
    gfrAccRA, bIsStatRA, qRA, wRA, ...
    dPelvis, dLFemur, dRFemur, dLTibia, dRTibia, uwb_mea, options)
    [nSamples, ~] = size(gfrAccMP);
    
    fOpt = struct('fs', 60);
    
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
    xhat_pos.RPV = zeros(4,4,nSamples);
    xhat_pos.RPV(1:3,1:3,:) = quat2rotm(qMP); 
    xhat_pos.RPV(4,4,:) = 1;
    xhat_pos.LSK = zeros(4,4,nSamples);
    xhat_pos.LSK(1:3,1:3,:) = quat2rotm(qLA); 
    xhat_pos.LSK(4,4,:) = 1;
    xhat_pos.RSK = zeros(4,4,nSamples);
    xhat_pos.RSK(1:3,1:3,:) = quat2rotm(qRA); 
    xhat_pos.RSK(4,4,:) = 1;

    pRPV = x0(01:03)'; vRPV = x0(04:06)';
    pLSK = x0(11:13)'; vLSK = x0(14:16)';
    pRSK = x0(21:23)'; vRSK = x0(24:26)';
    for n=1:nSamples
        xhat_pos.RPV(1:3,4,n) = (pRPV + vRPV*dt + gfrAccMP(n,:)*dt2)';
        xhat_pos.LSK(1:3,4,n) = (pLSK + vLSK*dt + gfrAccLA(n,:)*dt2)';
        xhat_pos.RSK(1:3,4,n) = (pRSK + vRSK*dt + gfrAccRA(n,:)*dt2)';
        
        vRPV = vRPV + gfrAccMP(n,:)*dt;
        vLSK = vLSK + gfrAccLA(n,:)*dt;
        vRSK = vRSK + gfrAccRA(n,:)*dt;
    end
        
    LTIBz = squeeze(xhat_pos.LSK(1:3,3,:))';
    RTIBz = squeeze(xhat_pos.RSK(1:3,3,:))';
    PELVy = squeeze(xhat_pos.RPV(1:3,2,:))';
    debug_dat.LFEO = squeeze(xhat_pos.LSK(1:3,4,:))' + dLTibia * LTIBz;
    debug_dat.RFEO = squeeze(xhat_pos.RSK(1:3,4,:))' + dRTibia * RTIBz;
    debug_dat.LFEP = squeeze(xhat_pos.RPV(1:3,4,:))' + dPelvis/2 * PELVy;
    debug_dat.RFEP = squeeze(xhat_pos.RPV(1:3,4,:))' - dPelvis/2 * PELVy;
    
    R_LFEM = zeros(3,3,nSamples);
    R_RFEM = zeros(3,3,nSamples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(xhat_pos.LSK(1:3,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(xhat_pos.RSK(1:3,2,:));
    RFEM_x = cross(RFEM_y, RFEM_z);
    R_RFEM(:,3,:) = RFEM_z ./ vecnorm(RFEM_z, 2, 1);
    R_RFEM(:,2,:) = RFEM_y ./ vecnorm(RFEM_y, 2, 1);
    R_RFEM(:,1,:) = RFEM_x ./ vecnorm(RFEM_x, 2, 1);
    debug_dat.qLTH = rotm2quat(R_LFEM);
    debug_dat.qRTH = rotm2quat(R_RFEM);
    
    xhat_pri = 0;
    xhat_pos.x = zeros(nSamples, 1);
end