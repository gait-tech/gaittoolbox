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

    fOpt = struct('fs', 60);
    
    xhat_pos.RPV
    xhat_pos.LSK
    xhat_pos.RSK
    xhat_pos.x
    
    
end