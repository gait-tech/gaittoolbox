% clc; clear all; close all;
addpath ./Visualize;


%% function [] = main()

% testfile = '50F0702.mat';
testfile = '50F0706.mat';
% testfile = '50F0107.mat';
% path containing vicon bone segment data 'Vicon Bone Segment Data from LK_MDR'
BASE_DIR = fullfile(pwd,'\BTK_EXPORT_MATCH\');

% path containing functions to transform vicon data to simulated imu
addpath(fullfile(BASE_DIR,'\ViconToIMU_functions\'));
fulltestfile = fullfile(BASE_DIR,testfile);
rawimport = importdata(fulltestfile);

%  ------------------------------------------------------------------------
%  Do not modify these variables
segmentMarkers={
'LFEO';'LFEA';'LFEP';'LFEL'; % the left femur
'LTIO';'LTIA';'LTIP';'LTIL'; % the left tibia
'PELO';'PELA';'PELP';'PELL'; % the pelvis
'RFEO';'RFEA';'RFEP';'RFEL'; % the right femur
'RTIO';'RTIA';'RTIP';'RTIL'; % the right tibia
'LFOO';'LFOA';'LFOP';'LFOL'; % the left foot
'RFOO';'RFOA';'RFOP';'RFOL'; % the right foot
}; 

sensorSegmentNames={'Left_Femur';'Left_Tibia';'Mid_Pelvis';
                'Right_Femur';'Right_Tibia';'Left_Foot';'Right_Foot';};

sensorSegmentOrigins={'LFEO';'LTIO';'PELO';'RFEO';'RTIO';'LFOO';'RFOO';};
%  ------------------------------------------------------------------------

fs = 100;
viconEvents = [];
tbl_vicon_events = [];
if length( fieldnames(rawimport) ) == 2
    markers = rawimport.viconmarkers;
    viconEvents = rawimport.eventsIndex; 
    name_events = fieldnames(viconEvents);
    N_EVENTS = length(name_events);
    boolViconEvents = false(length(markers.LANK(:,1)),4);
    tbl_vicon_events = ...
        array2table(boolViconEvents,...
        'VariableNames',name_events);
    for e = 1:N_EVENTS
        c_event = tbl_vicon_events.(name_events{e});
        c_event(viconEvents.(name_events{e})) = true;
        tbl_vicon_events.(name_events{e}) = c_event;
    end
else
    markers = rawimport;
end

bMarkersMissing = ~isfield(markers,segmentMarkers);
if sum(bMarkersMissing) > 0, 
    fprintf('Markers missing in %s\n',testfile);
    return; % break
end
% crop all rows missing markers
tbl_markers = struct2table(markers);
arr_markers = table2array(tbl_markers)./1000;
[NROWS,~]   = find(arr_markers==0);
unique_rows = unique(NROWS);
% remove these rows from the table
tbl_markers(unique_rows,:) = [];
tbl_vicon_events(unique_rows,:) = [];
arr_markers(unique_rows,:) = [];
[~,NCOLS] = size(tbl_markers);

% iterate over columns and convert from mm to m
tbl_markers{:,:} = tbl_markers{:,:}./1000;
% add padding to give kf some time to converge
% tbl_markers = [tbl_markers(ones(299,1),:); tbl_markers];

% specify noise level in the accelerometer (m/s^2)
sigma_acc = 0.5;
% specify noise level in the gyro (rad/s)
sigma_gyr = 0.05;
% specify noise level in the gyro (microT)
sigma_mag = 0.5;

% specify reference measurement of gravity
gRef = 9.796720; %http://iopscience.iop.org/article/10.1088/0026-1394/9/2/001/meta

% define the rotation matrix which converts all of the measurements in the
% vicon frame to the 'psuedo' global frame of reference
qViconGfr = [1 0 0 0];

% --- low-pass filter the vicon trajectories to reduce marker flicker
fc    = 10;
[b,i] = butter(2,fc/(fs/2),'low');
N_MARKERS = length(segmentMarkers);
for m = 1:N_MARKERS
    tbl_markers.(segmentMarkers{m}) = ...
        filtfilt(b,i,tbl_markers.(segmentMarkers{m}));
end

mRef = [sqrt(24157^2+5381^2), 0, -51428]./1000; % microT -> 1 micro = 1000 nano
% Australian Geomagnetic Reference Field Computation
% Requested: Latitude -33o 55' 8", Longitude 151o 13' 52", Elevation 0 km, Date 2017/07/6 
% Calculated: Latitude -33.9189o, Longitude +151.2311o, Elevation 0.00 km, Epoch 2017.5096
% D = 12.558 deg | H = 24749 nT | I = -64.301 deg 
% X = 24157 nT   | Y = 5381 nT  | Z = -51428 nT |

[N_SEGMENTS,~] = size(sensorSegmentNames);
imuSegment  = struct;
quatSegment = struct;
gfrSegment  = struct;
%% ---------------------------------------------------------------------
%  --- Simulated MIMU Measurements from body segments

% --- Optional, add local magnetic disturbance
% random orientation
u1 = rand(1); u2 = rand(1); u3 = rand(1);
% http://planning.cs.uiuc.edu/node198.html
randomQuaternion = ... 
quatnormalize([sqrt(1-u1)*sin(2*pi*u2),...
               sqrt(1-u1)*cos(2*pi*u2),...
               sqrt(1-u1)*sin(2*pi*u3),...
               sqrt(1-u1)*cos(2*pi*u3)]);
           
dipoleOrigin = [0; 0; 0];
magneticFieldStrength = 5;% microTesla
dipoleObj = dipole( quaternion2(randomQuaternion') , dipoleOrigin, ...
    magneticFieldStrength);
for r = 1:N_SEGMENTS
    oPos = tbl_markers.(segmentMarkers{ (r-1)*4 + 1} );
    aPos = tbl_markers.(segmentMarkers{ (r-1)*4 + 2} );
    pPos = tbl_markers.(segmentMarkers{ (r-1)*4 + 3} );
    lPos = tbl_markers.(segmentMarkers{ (r-1)*4 + 4} );

    [quat,imu,gfr] = viconGenIMU('fs',fs,...
        'origin',oPos,'anterior',aPos,'proximal',pPos,'lateral',lPos,...
        'sigma_acc',0,'sigma_mag',sigma_mag,'sigma_gyr',sigma_gyr,...
        'qViconGfr',qViconGfr,'gRef',gRef,'mRef',mRef,'dipoleOrigin',dipoleObj,...
        'bOutputTable',false);    

    quatSegment.(sensorSegmentNames{r}) = quat;
    imuSegment.(sensorSegmentNames{r}) = imu;
    gfrSegment.(sensorSegmentNames{r}) = gfr;
end

%% -----------------------------------------------------------------------
%  Simulate uwb measurement by generating pairwise combinations, using the
%  origin of each bone segment as the root point
uwb_mea = struct;
sensorSegmentPairs = nchoosek(1:N_SEGMENTS,2);
for u = 1:length(sensorSegmentPairs(:,1))
    p1 = sensorSegmentPairs(u,1);
    o1 = tbl_markers.(segmentMarkers{ (p1-1)*4 + 1} );
    s1 = lower(sensorSegmentNames{p1}); 
    
    p2 = sensorSegmentPairs(u,2);
    o2 = tbl_markers.(segmentMarkers{ (p2-1)*4 + 1} );
    s2 = lower(sensorSegmentNames{p2}); 
    
    euclid_dist = vecnormalize(o2-o1);
    uwb_fieldname = strcat(s1,'_',s2);
    % add measurement as field in struct
    uwb_mea.(uwb_fieldname) = euclid_dist;
end

%% -----------------------------------------------------------------------
%  compare the magnetometer measurements with and without magnetic
%  interference
% figure('name','Simulated Magnetometer Measurements GFR');
% hs(1)=subplot(3,1,1);
% plot([gfrSegment.Left_Foot.mag])
% title('Local Magnetic Field');
% hs(2)=subplot(3,1,2);
% plot([gfrSegment.Left_Foot.mag_dis])
% title('Local Magnetic Field Disturbance');
% hs(3)=subplot(3,1,3);
% plot([gfrSegment.Left_Foot.mag + gfrSegment.Left_Foot.mag_dis]);
% title('Local Magnetic Field + Disturbance');
% linkaxes(hs,'x');

% figure('name','Simulated Magnetometer Measurements IMU');
% hs(1)=subplot(3,1,1);
% plot([imuSegment.Left_Foot.mag])
% title('Local Magnetic Field');
% hs(2)=subplot(3,1,2);
% plot([imuSegment.Left_Foot.mag_dis])
% title('Local Magnetic Field Disturbance');
% hs(3)=subplot(3,1,3);
% plot([imuSegment.Left_Foot.mag + imuSegment.Left_Foot.mag_dis]);
% title('Local Magnetic Field + Disturbance');
% linkaxes(hs,'x');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
%% Start reconstructing the lower body
% --- estimate orientation of each sensor and express acceleration in GFR

% local variable assignment for readability
imuPelvis = imuSegment.Mid_Pelvis;
imuLankle = imuSegment.Left_Tibia;
imuRankle = imuSegment.Right_Tibia;

gfrPelvis = gfrSegment.Mid_Pelvis;
gfrLankle = gfrSegment.Left_Tibia;
gfrRankle = gfrSegment.Right_Tibia;

qPelvis = quatSegment.Mid_Pelvis;
qLankle = quatSegment.Left_Tibia;
qRankle = quatSegment.Right_Tibia;

N_TIME_PTS = length(qPelvis(:,1));
muAcc = 0.003;
muMag = 0.001;

qPelvisEst = wrapper_MIMU_CAHRS_ArcTan('qInit',qPelvis(1,:),...
    'fs',fs,'Acc',imuPelvis.acc,'Gyr',imuPelvis.gyr,....
    'Mag',imuPelvis.mag,...
    'static_mu_acc',0.5,'dynamic_mu_acc',muAcc,...
    'static_mu_mag',0.5,'dynamic_mu_mag',muMag,...
    'bDynamicMu',true(N_TIME_PTS,1));
[rmea,pmea,ymea]=quaternion2nautical(qPelvis,'');
[rdeg,pdeg,ydeg]=quaternion2nautical(qPelvisEst,'');
% %%
[s_xhatMP,s_yhatMP,s_zhatMP] = quaternion2rotationmatrix(qPelvis,'column');

sigma_A = 0.01;
sigma_G = 0.01;
sigma_M = 0.1;
% [ xPri,xPos,Ppri,Ppos,sAcc ] = ...
%     parkKF( fs, imuPelvis.acc, sigma_A, imuPelvis.gyr, sigma_G,...
%     0.1 ,s_zhatMP(1,:) );

h_ref = gfrPelvis.mag(1,:);
g_ref = gfrPelvis.acc(1,:) - gfrPelvis.bod_acc(1,:);
% rotateVecsByQuats(imuPelvis.acc(1,:)-imuPelvis.body_acc(1,:),qPelvis(1,:));

[~,dip_ref_rad] = angle_two_vectors(g_ref,  h_ref);
if dip_ref_rad > pi/2
    dip_ref_rad = dip_ref_rad - pi/2;
end
dip_ref_deg = rad2deg(dip_ref_rad);


% [xk_up_pri,xk_up_pos,Pk_up_pri,Pk_up_pos,sAcc,...
%     xk_yaw_pri,xk_yaw_pos,Pk_yaw_pri,Pk_yaw_pos,qPelvisKF_est] ...
%     = cascaded_AHRS( fs, imuPelvis.acc, sigma_A, imuPelvis.gyr, sigma_G,...
%     0.3 ,s_zhatMP(1,:), [imuPelvis.mag + imuPelvis.mag_dis],sigma_M,s_xhatMP(1,:), h_ref, dip_ref_rad);
% roll  = atan( xk_up_pos(:,2)./ xk_up_pos(:,3)); % gamma 
% rolldeg = rad2deg(roll);
% pitch = atan(-xk_up_pos(:,1)./(xk_up_pos(:,2)./sin(roll))); % beta
% pitchdeg = rad2deg(pitch);
% yaw = atan2( (-cos(roll).*xk_yaw_pos(:,2) + sin(roll).*xk_yaw_pos(:,3)),...
%     (xk_yaw_pos(:,1)./cos(pitch)) );
% yawdeg = rad2deg(yaw);

% [rolldeg,pitchdeg,yawdeg] = quaternion2nautical(qPelvisKF_est,'');

% updateFigureContents('Attitude Comparison');
% subplot(3,1,1);hold on;
% plot(rmea,'-k','LineWidth',2);
% plot(rdeg,':xb');
% plot(rolldeg,':+r');
% subplot(3,1,2);hold on;
% plot(pmea,'-k','LineWidth',2);
% plot(pdeg,':xb');
% plot(pitchdeg,':+r');
% subplot(3,1,3);hold on;
% plot(ymea,'-k','LineWidth',2);
% plot(ydeg,':xb');
% plot(yawdeg,':+r');
% 
% updateFigureContents('pwelch accelerometer magnetometer');
% subplot(2,1,1);
% pwelch(vecnormalize(imuPelvis.acc),[],[],[],fs,'onesided');
% subplot(2,1,2);
% pwelch(vecnormalize(imuLankle.mag+imuLankle.mag_dis),[],[],[],fs,'onesided');
% break
bNeg = qPelvisEst(:,1) < 0;
qPelvisEst(bNeg,:) = -qPelvisEst(bNeg,:);
% gfrMagPelvisTru=rotateVecsByQuats(imuPelvis.mag+imuPelvis.mag_dis,qPelvisEst);
gfrMagPelvisEst=rotateVecsByQuats(imuPelvis.mag+imuPelvis.mag_dis,qPelvisEst);
% gfrMagPelvis2=quatmultiply(...
%     quatmultiply(quatnormalize(qPelvisEst),...
%     [zeros(N_TIME_PTS,1), imuPelvis.mag]), quatinv(quatnormalize(qPelvisEst)));

% Plot pelvis euler angles and quaternions
% figure('name','euler angles');hold on;
% plot(quaternion2nautical(qPelvis),':');
% plot(quaternion2nautical(qPelvisEst),'*');
% 
% figure('name','Quaternions');
% subplot(4,1,1);hold on;
% plot((qPelvis(:,1)),':');
% plot((qPelvisEst(:,1)),'*');
% 
% subplot(4,1,2);hold on;
% plot((qPelvis(:,2)),':');
% plot((qPelvisEst(:,2)),'*');
% 
% subplot(4,1,3);hold on;
% plot((qPelvis(:,3)),':');
% plot((qPelvisEst(:,3)),'*');
% 
% subplot(4,1,4);hold on;
% plot((qPelvis(:,4)),':');
% plot((qPelvisEst(:,4)),'*');


gfr_acc_MP = rotateVecsByQuats(imuPelvis.acc,qPelvisEst );
gfrPelvisEstGyr = rotateVecsByQuats(imuPelvis.gyr,qPelvisEst);

qLankleEst = wrapper_MIMU_CAHRS_ArcTan('qInit',qLankle(1,:),...
    'fs',fs,'Acc',imuLankle.acc,'Gyr',imuLankle.gyr,'Mag',imuLankle.mag,...
    'static_mu_acc',0.5,'dynamic_mu_acc',muAcc,...
    'static_mu_mag',0.5,'dynamic_mu_mag',muMag,...
    'bDynamicMu',true(N_TIME_PTS,1));

qLankleEstMagDist = wrapper_MIMU_CAHRS_ArcTan('qInit',qLankle(1,:),...
    'fs',fs,'Acc',imuLankle.acc,'Gyr',imuLankle.gyr,...
    'Mag',imuLankle.mag,...
    'static_mu_acc',0.5,'dynamic_mu_acc',muAcc,...
    'static_mu_mag',0.5,'dynamic_mu_mag',muMag,...
    'bDynamicMu',true(N_TIME_PTS,1));

gfrLankleEstMag = rotateVecsByQuats( [imuLankle.mag+imuLankle.mag_dis],...
    qLankleEst );

mag_perb_detect = (gfrLankleEstMag-gfrLankle.mag);
gyr_norm = sqrt(sum(imuLankle.gyr.^2,2));

gfr_acc_LA = rotateVecsByQuats( imuLankle.acc,qLankleEst );
gfrLankleEstGyr = rotateVecsByQuats( imuLankle.gyr,qLankleEst );

qRankleEst = wrapper_MIMU_CAHRS_ArcTan('qInit',qRankle(1,:),...
    'fs',fs,'Acc',imuRankle.acc,'Gyr',imuRankle.gyr,'Mag',imuRankle.mag,...
    'static_mu_acc',0.5,'dynamic_mu_acc',muAcc,...
    'static_mu_mag',0.5,'dynamic_mu_mag',muMag,...
    'bDynamicMu',true(N_TIME_PTS,1));
gfr_acc_RA = rotateVecsByQuats(imuRankle.acc,qRankleEst);
gfrRankleEstGyr = rotateVecsByQuats(imuRankle.gyr,qRankleEst);

% subtract gravity
gfr_acc_MP(:,3) = gfr_acc_MP(:,3)-gRef;
gfr_acc_LA(:,3) = gfr_acc_LA(:,3)-gRef;
gfr_acc_RA(:,3) = gfr_acc_RA(:,3)-gRef;

% -------------------------------------------------------------------------
%% detect stationary periods thresholding on the norm of accelerometer
WIN_SECS = 0.25;
VAR_WIN  = floor(fs*WIN_SECS); % NUM_SAMPLES
ACC_VAR_THRESH = 1;

movVarAcc_pelvis = movingvar(sqrt( sum(imuPelvis.acc.^2,2)) ,VAR_WIN);
bIsStatMP = movVarAcc_pelvis < 0;
accMagMP = vecnormalize(imuPelvis.acc);
gyrMagMP = vecnormalize(imuPelvis.gyr);

movVarAcc_lankle = movingvar(sqrt( sum(imuLankle.acc.^2,2)) ,VAR_WIN);
bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
accMagLA = vecnormalize(imuLankle.acc);
gyrMagLA = vecnormalize(imuLankle.gyr);

movVarAcc_rankle = movingvar(sqrt( sum(imuRankle.acc.^2,2)) ,VAR_WIN);
bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
accMagRA = vecnormalize(imuRankle.acc);
gyrMagRA = vecnormalize(imuRankle.gyr);

% Plot Stationary Detector
% figure('name','Stationary Detector');
% subplot(2,1,1);hold on;hp=[];lgstr={};
% hp(end+1)=plot(accMagMP,'k');lgstr{end+1}='mid-pelvis';
% hp(end+1)=plot(accMagLA,'b');lgstr{end+1}='left ankle';
% hp(end+1)=plot(accMagRA,'r');lgstr{end+1}='right ankle';
% legend(hp,lgstr);
% 
% subplot(2,1,2);hold on;hp=[];lgstr={};
% hp(end+1)=plot(movVarAcc_pelvis,'k');lgstr{end+1}='mid-pelvis';
% hp(end+1)=plot(movVarAcc_lankle,'b');lgstr{end+1}='left ankle';
% hp(end+1)=plot(movVarAcc_rankle,'r');lgstr{end+1}='right ankle';
% legend(hp,lgstr);

% define boolean vectors for x,y,z components of zero velocity update
bIsStationaryMP = repmat(bIsStatMP,1,3);
bIsStationaryLA = repmat(bIsStatLA,1,3);
bIsStationaryRA = repmat(bIsStatRA,1,3);

% assume start position of each sensor is known
x0_pos_MP = tbl_markers.PELO(1,:);
x0_pos_LA = tbl_markers.LTIO(1,:);
x0_pos_RA = tbl_markers.RTIO(1,:);
x0_vel_MP = gfrPelvis.vel(1,:);
x0_vel_LA = gfrLankle.vel(1,:);
x0_vel_RA = gfrRankle.vel(1,:);
% -------------------------------------------------------------------------
%% ------- Reconstruct mid pelvis and ankle joint centre trajectories -----
% Stephen add your code within this wrapper function

[xMP,PxMP,yMP,PyMP,zMP,PzMP] = ...
    wrapper_PosVel_3KFs( fs,gfr_acc_MP,...
    [x0_pos_MP, x0_vel_MP]',1,bIsStationaryMP );

[xLA,PxLA,yLA,PyLA,zLA,PzLA] = ...
    wrapper_PosVel_3KFs( fs,gfr_acc_LA,...
    [x0_pos_LA, x0_vel_LA]',1,bIsStationaryLA);

[xRA,PxRA,yRA,PyRA,zRA,PzRA] = ...
    wrapper_PosVel_3KFs( fs,gfr_acc_RA,...
    [x0_pos_RA, x0_vel_RA]',1,bIsStationaryRA);

% [xMP,PxMP,yMP,PyMP,zMP,PzMP] = ...
%     wrapper_PosVel_3KFs( fs,gfrSegment.Mid_Pelvis.bod_acc,...
%     [x0_pos_MP, x0_vel_MP]',1,bIsStationaryMP );
% 
% [xLA,PxLA,yLA,PyLA,zLA,PzLA] = ...
%     wrapper_PosVel_3KFs( fs,gfrSegment.Left_Tibia.bod_acc,...
%     [x0_pos_LA, x0_vel_LA]',1,bIsStationaryLA);
% 
% [xRA,PxRA,yRA,PyRA,zRA,PzRA] = ...
%     wrapper_PosVel_3KFs( fs,gfrSegment.Right_Tibia.bod_acc,...
%     [x0_pos_RA, x0_vel_RA]',1,bIsStationaryRA);
% -------------------------------------------------------------------------
% Alternative approach with UWB measurements

% [ x_pri,x_pos ] = kf_3_kmus(fs, sigma_acc, ...
%     x0_pos_MP, x0_vel_MP, gfrSegment.Mid_Pelvis.bod_acc, bIsStationaryMP,...
%     x0_pos_LA, x0_vel_LA, gfrSegment.Left_Tibia.bod_acc, bIsStationaryLA,...
%     x0_pos_RA, x0_vel_RA, gfrSegment.Right_Tibia.bod_acc, bIsStationaryRA, uwb_mea);

[ x_pri, x_pos ] = kf_3_kmus(fs, sigma_acc, ...
    x0_pos_MP, x0_vel_MP, gfr_acc_MP, bIsStatMP,...
    x0_pos_LA, x0_vel_LA, gfr_acc_LA, bIsStatLA,...
    x0_pos_RA, x0_vel_RA, gfr_acc_RA, bIsStatRA, uwb_mea);

d_pelvis = norm(tbl_markers.RFEP(1,:)-tbl_markers.LFEP(1,:));
d_rfemur = norm(tbl_markers.RFEP(1,:)-tbl_markers.RFEO(1,:));
d_lfemur = norm(tbl_markers.LFEP(1,:)-tbl_markers.LFEO(1,:));
d_rtibia = norm(tbl_markers.RFEO(1,:)-tbl_markers.RTIO(1,:));
d_ltibia = norm(tbl_markers.LFEO(1,:)-tbl_markers.LTIO(1,:));

[ x_pri_v2, x_pos_v2, t_dat_v2 ] = kf_3_kmus_v2(fs, ...
    sigma_acc, sigma_acc, sigma_acc, false, ...
    x0_pos_MP, x0_vel_MP, gfr_acc_MP, bIsStatMP, qPelvisEst, ...
    x0_pos_LA, x0_vel_LA, gfr_acc_LA, bIsStatLA, qLankleEst, ...
    x0_pos_RA, x0_vel_RA, gfr_acc_RA, bIsStatRA, qRankleEst, ...
    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, uwb_mea, ...
    true, false, true, false, false);

% ------------------------------------------------------------------------
% Visualise Trajectories for sanity check
% figure('name','Animation Vicon Bone Segments');
% xlabel('x - Forward (m)');ylabel('y - East (m)');zlabel('z - Vertical (m)');
% hold on;grid on;axis('equal');view([-51 12]);%view([0 0]);

PELO = tbl_markers.PELO; PELP = tbl_markers.PELP;

LFEO = tbl_markers.LFEO; LFEP = tbl_markers.LFEP;
LTIO = tbl_markers.LTIO; LTIP = tbl_markers.LTIP;
LFOO = tbl_markers.LFOO; LFOP = tbl_markers.LFOP;

RFEO = tbl_markers.RFEO; RFEP = tbl_markers.RFEP;
RTIO = tbl_markers.RTIO; RTIP = tbl_markers.RTIP;
RFOO = tbl_markers.RFOO; RFOP = tbl_markers.RFOP;

PMid = [mean([RFEP(:,1), LFEP(:,1)],2),...
        mean([RFEP(:,2), LFEP(:,2)],2),...
        mean([RFEP(:,3), LFEP(:,3)],2)];

% % modify this if you want to focus on a smaller time period of the
% % trajectory
% allIdx = 1:N_TIME_PTS; 
% arr_markers = [table2array(tbl_markers), ...
%     x_pos(:,[1:3,7:9,13:15]),xMP,yMP,zMP, xLA,yLA,zLA, xRA,yRA,zRA];
% 
% vicX = reshape(arr_markers(allIdx,1:3:end),[],1); XLIM = [min(vicX),max(vicX)];
% vicY = reshape(arr_markers(allIdx,2:3:end),[],1); YLIM = [min(vicY),max(vicY)];
% vicZ = reshape(arr_markers(allIdx,3:3:end),[],1); ZLIM = [min(vicZ),max(vicZ)];
% xlim(XLIM);ylim(YLIM);zlim(ZLIM);
% 
% % KF estimates
% p1 = plot3(xMP(allIdx,1),yMP(allIdx,1),zMP(allIdx,1),'.k');
% plot3(xLA(allIdx,1),yLA(allIdx,1),zLA(allIdx,1),'.b');
% plot3(xRA(allIdx,1),yRA(allIdx,1),zRA(allIdx,1),'.r');
% 
% % % KF with UWB
% % p2 = plot3(x_pos(allIdx,1), x_pos(allIdx,2), x_pos(allIdx,3),'ok');
% % plot3(x_pos(allIdx,7), x_pos(allIdx,8), x_pos(allIdx,9),'oc');
% % plot3(x_pos(allIdx,13),x_pos(allIdx,14),x_pos(allIdx,15),'om');
% 
% % KF v2
% p3 = plot3(x_pos_v2(allIdx,1), x_pos_v2(allIdx,2), x_pos_v2(allIdx,3),'+k');
% plot3(x_pos_v2(allIdx,7), x_pos_v2(allIdx,8), x_pos_v2(allIdx,9),'+c');
% plot3(x_pos_v2(allIdx,13),x_pos_v2(allIdx,14),x_pos_v2(allIdx,15),'+m');
% 
% % TRUE positions
% h_mp_mea = plot3(PMid(allIdx,1),PMid(allIdx,2),PMid(allIdx,3),'-k');
% h_la_mea = plot3(LTIO(allIdx,1),LTIO(allIdx,2),LTIO(allIdx,3),'-b');
% h_ra_mea = plot3(RTIO(allIdx,1),RTIO(allIdx,2),RTIO(allIdx,3),'-r');
% 
% % updateFigureContents('Animation');
% 
% % %%
% % updateFigureContents('Animation Vicon Bone Segments');
% % xlabel('x - Forward (m)');ylabel('y - East (m)');zlabel('z - Vertical (m)');
% % hold on;grid on;axis('equal');view([-51 12]);%view([0 0]);
% % xlim(XLIM);ylim(YLIM);zlim(ZLIM);
% % 
% % fctr=1;
% % 
% % % for i=N_TIME_PTS
% %     
% % cla;
% i=N_TIME_PTS;
% % -------- Pelvis --------
% plot3(PELO(i,1),PELO(i,2),PELO(i,3),'ok','MarkerFaceColor','k');
% plot3(PELP(i,1),PELP(i,2),PELP(i,3),'ok','MarkerFaceColor','none');
% line([PELO(i,1), PELP(i,1)],...
%      [PELO(i,2), PELP(i,2)],...
%      [PELO(i,3), PELP(i,3)],...
% 'Color','k','LineWidth',2);
% line([LFEP(i,1), RFEP(i,1)],...
%      [LFEP(i,2), RFEP(i,2)],...
%      [LFEP(i,3), RFEP(i,3)],...
% 'Color','k','LineWidth',3);
% 
% % Left Femur
% plot3(LFEO(i,1),LFEO(i,2),LFEO(i,3),'ob','MarkerFaceColor','b');
% plot3(LFEP(i,1),LFEP(i,2),LFEP(i,3),'ob','MarkerFaceColor','none');
% line([LFEO(i,1),LFEP(i,1)],[LFEO(i,2), LFEP(i,2)],[LFEO(i,3), LFEP(i,3)],...
% 'Color','c','LineWidth',2,'LineStyle','--');
% 
% % Left Tibia
% plot3(LTIO(i,1),LTIO(i,2),LTIO(i,3),'ob','MarkerFaceColor','b');
% plot3(LTIP(i,1),LTIP(i,2),LTIP(i,3),'ob','MarkerFaceColor','none');
% line([LTIO(i,1),LTIP(i,1)],[LTIO(i,2), LTIP(i,2)],[LTIO(i,3), LTIP(i,3)],...
% 'Color','c','LineWidth',2);
% 
% % Left Foot
% plot3(LFOO(i,1),LFOO(i,2),LFOO(i,3),'ob','MarkerFaceColor','b');
% plot3(LFOP(i,1),LFOP(i,2),LFOP(i,3),'ob','MarkerFaceColor','none');
% line([LFOO(i,1),LFOP(i,1)],[LFOO(i,2), LFOP(i,2)],[LFOO(i,3), LFOP(i,3)],...
% 'Color','c','LineWidth',2,'LineStyle','--');
% 
% % Right Femur
% plot3(RFEO(i,1),RFEO(i,2),RFEO(i,3),'or','MarkerFaceColor','r');
% plot3(RFEP(i,1),RFEP(i,2),RFEP(i,3),'or','MarkerFaceColor','none');
% line([RFEO(i,1), RFEP(i,1)],[RFEO(i,2), RFEP(i,2)],[RFEO(i,3), RFEP(i,3)],...
% 'Color','m','LineWidth',2,'LineStyle','--');
% 
% % Right Tibia
% plot3(RTIO(i,1),RTIO(i,2),RTIO(i,3),'or','MarkerFaceColor','r');
% plot3(RTIP(i,1),RTIP(i,2),RTIP(i,3),'or','MarkerFaceColor','none');
% line([RTIO(i,1), RTIP(i,1)],[RTIO(i,2), RTIP(i,2)],[RTIO(i,3), RTIP(i,3)],...
% 'Color','m','LineWidth',2);
% 
% % Right Foot
% plot3(RFOO(i,1),RFOO(i,2),RFOO(i,3),'or','MarkerFaceColor','r');
% plot3(RFOP(i,1),RFOP(i,2),RFOP(i,3),'or','MarkerFaceColor','none');
% line([RFOO(i,1),RFOP(i,1)],[RFOO(i,2), RFOP(i,2)],[RFOO(i,3), RFOP(i,3)],...
% 'Color','m','LineWidth',2,'LineStyle','--');
% % 
% % title(['Time:= ', num2str(i/fs)]);
% % % drawnow;
% % % F(fctr) = getframe(gcf);
% % fctr = fctr+1;
% % 
% % pause(0.5);
% % end
% hand.fill2=fill([-3,3,3,-3],[-4,-4,2,2],[0,0,0,0]);
% set(hand.fill2,'facealpha',0.51);
% set(hand.fill2,'edgealpha',0.1);
% 
% legend('off')
% legend([p1, p3, h_mp_mea], 'ekf', 'ekf v2', 'true pos');

% video = VideoWriter('ViconWalk.avi','Motion JPEG AVI');
% open(video);
% writeVideo(video,F);
% close(video);

% ------------------------------------------------------------------------
% Result Validations / Debugging
idx = 1:length(qPelvisEst(:,1));

estBody = Body('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
               'lnSymbol', '--', 'ptSymbol', 'o', ...
               'SACR', x_pos_v2(idx,1:3), 'LFEP', t_dat_v2.LFEP(idx,:), ...
               'LFEO', t_dat_v2.LFEO(idx,:), 'LTIO', x_pos_v2(idx,7:9), ...
               'RFEP', t_dat_v2.RFEP(idx,:), ...
               'RFEO', t_dat_v2.RFEO(idx,:), ...
               'RTIO', x_pos_v2(idx,13:15), ...
               'qPelvis', qPelvisEst(idx,:), ...
               'qRFemur', t_dat_v2.qRFemur(idx,:), ...
               'qLFemur', t_dat_v2.qLFemur(idx,:), ...
               'qRTibia', qRankleEst(idx,:), ...
               'qLTibia', qLankleEst(idx,:));
           
actBody = Body('name', 'act', 'posUnit', 'm', 'oriUnit', 'deg', ...
               'lnSymbol', '-', 'ptSymbol', '.', ...
               'SACR', PMid(idx,:), 'LFEP', LFEP(idx,:), ...
               'LFEO', LFEO(idx,:), 'LTIO', LTIO(idx,:), ...
               'RFEP', RFEP(idx,:), 'RFEO', RFEO(idx,:), ...
               'RTIO', RTIO(idx,:), ...
               'qPelvis', quatSegment.Mid_Pelvis(idx,:), ...
               'qRFemur', quatSegment.Right_Femur(idx,:), ...
               'qLFemur', quatSegment.Left_Femur(idx,:), ...
               'qRTibia', quatSegment.Right_Tibia(idx,:), ...
               'qLTibia', quatSegment.Left_Tibia(idx,:));

actMPVel = [0 0 0; actBody.SACR(2:end, :)-actBody.SACR(1:end-1, :)]*fs;
actLAVel = [0 0 0; actBody.LTIO(2:end, :)-actBody.LTIO(1:end-1, :)]*fs;
actRAVel = [0 0 0; actBody.RTIO(2:end, :)-actBody.RTIO(1:end-1, :)]*fs;
actState = [PMid actMPVel LTIO actLAVel RTIO actRAVel];

% % Position error
% updateFigureContents('Position (1)');
% plotPosition({estBody, actBody}, {'SACR'});
% updateFigureContents('Position (2)');
% plotPosition({estBody, actBody}, {'LTIO', 'RTIO'});
% updateFigureContents('Position (3)');
% plotPosition({estBody, actBody}, {'LFEO', 'RFEO'});
% 
% updateFigureContents('Position Error');
% plotPositionDiff(actBody, estBody, {'SACR', 'LTIO', 'RTIO'});
% 
% % Orientation Error
% updateFigureContents('Orientation (1)');
% plotOrientation({estBody, actBody}, {'qPelvis', 'qRTibia', 'qLTibia'});
% updateFigureContents('Orientation (2)');
% plotOrientation({estBody, actBody}, {'qRFemur', 'qLFemur'});
% 
% updateFigureContents('Orientation Error');
% plotOrientationDiff(actBody, estBody, {'qRFemur', 'qLFemur'});
% 
% % Segment Length Error
% updateFigureContents('Segment Length Error');
% plotLowerBodySegmentLengthError(estBody, d_pelvis, d_lfemur, d_rfemur, ...
%     d_ltibia, d_rtibia)
% 
% % Position and Velocity
% updateFigureContents('State Comparison');
% plotPosition({estBody, actBody}, {'SACR', 'LTIO', 'RTIO'});
% 
% updateFigureContents('Pelvis Velocity Comparison'); hold on;
% for i=1:3
%     subplot(3,1,i); hold on;
%     title(strcat('Pelvis ', 'w'+i));
%     plot(idx, x_pos_v2(idx, 3+i), ...
%          strcat(estBody.xyzColor{1}, estBody.lnSymbol));
%     plot(idx, actMPVel(idx,i), ...
%          strcat(actBody.xyzColor{1}, actBody.lnSymbol));
% end
% 
% updateFigureContents('LTIO Velocity Comparison'); hold on;
% for i=1:3
%     subplot(3,1,i); hold on;
%     title(strcat('LTIO ', 'w'+i));
%     plot(idx, x_pos_v2(idx, 9+i), ...
%          strcat(estBody.xyzColor{1}, estBody.lnSymbol));
%     plot(idx, actLAVel(idx,i), ...
%          strcat(actBody.xyzColor{1}, actBody.lnSymbol));
% end
% 
% updateFigureContents('RTIO Velocity Comparison'); hold on;
% for i=1:3
%     subplot(3,1,i); hold on;
%     title(strcat('RTIO ', 'w'+i));
%     plot(idx, x_pos_v2(idx, 15+i), ...
%          strcat(estBody.xyzColor{1}, estBody.lnSymbol));
%     plot(idx, actRAVel(idx,i), ...
%          strcat(actBody.xyzColor{1}, actBody.lnSymbol));
% end

% % state progress
updateFigureContents('State Progress');
plotStateComparison(t_dat_v2, actState, 10);

% Why are there instances of "big changes" in the velocity of the MP, LA,
% and RA?

% % Snapshots
% updateFigureContents('Lower Body'); grid on;
% xlabel('x'); ylabel('y'); zlabel('z');
% for i=1:20:60
%     plotLowerBody(estBody, i);
% end
% 

% % Animation
% updateFigureContents('Animation');
% estBodyLimits = [estBody.xlim() estBody.ylim() estBody.zlim()];
% grid
% for i=1:613
%     clf;
%     xlim(estBodyLimits(1:2)); 
%     ylim(estBodyLimits(3:4)); 
%     zlim(estBodyLimits(5:6));  
%     view(0, 0);
%     plotLowerBody(estBody, i);
%     pause(1/10);
% end