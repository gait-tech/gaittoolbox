vicon = ViconNexus();

subj = 'S03'; 
[x1, y1, z1, e1] = vicon.GetTrajectory(subj, 'RASI');
[x2, y2, z2, e2] = vicon.GetTrajectory(subj, 'RPSI');

% firstIdx = max(find(e2, 1), 2);
% v = [x2(firstIdx) - x1(firstIdx), y2(firstIdx) - y1(firstIdx), ...
%      z2(firstIdx) - z1(firstIdx)];
% firstIdx = 1132;
v = [x2(end) - x1(end), y2(end) - y1(end), ...
     z2(end) - z1(end)];
x2(1:firstIdx-1) = x1(1:firstIdx-1) + v(1);
y2(1:firstIdx-1) = y1(1:firstIdx-1) + v(2);
z2(1:firstIdx-1) = z1(1:firstIdx-1) + v(3);
e2(1:firstIdx-1) = e2(firstIdx);

vicon.SetTrajectory(subj, 'RPSI', x2, y2, z2, e2);

%% rigid body fill (backwards)
subj = 'S10';
[x, y, z, e] = vicon.GetTrajectory(subj, 'LASI');
% forward
baseIdx = max(find(e, 1), 2);
targetIdx = 1:baseIdx-1;

% backward
baseIdx = find(e, 1, 'last');
targetIdx = baseIdx+1:length(e);

baseLASI = getPoints(subj, 'LASI', baseIdx);
baseRASI = getPoints(subj, 'RASI', baseIdx);
baseLPSI = getPoints(subj, 'LPSI', baseIdx);
baseRPSI = getPoints(subj, 'RPSI', baseIdx);

targetRASI = getPoints(subj, 'RASI', targetIdx);
targetLPSI = getPoints(subj, 'LPSI', targetIdx);
targetRPSI = getPoints(subj, 'RPSI', targetIdx);
targetLASI = rigidBodyFill(baseLPSI, baseRPSI, baseRASI, baseLASI, ...
                           targetLPSI, targetRPSI, targetRASI);
LASI = getPoints(subj, 'LASI');
LASI(targetIdx, :) = targetLASI;
setPoints(subj, 'LASI', LASI);

%% rigid body fill (backwards, RASI)
subj = 'S04';
[x, y, z, e] = vicon.GetTrajectory(subj, 'RASI');
% forward
baseIdx = max(find(e, 1), 2);
targetIdx = 1:baseIdx-1;

% backward
baseIdx = find(e, 1, 'last');
targetIdx = baseIdx+1:length(e);

baseLASI = getPoints(subj, 'LASI', baseIdx);
baseRASI = getPoints(subj, 'RASI', baseIdx);
baseLPSI = getPoints(subj, 'LPSI', baseIdx);
baseRPSI = getPoints(subj, 'RPSI', baseIdx);

targetLASI = getPoints(subj, 'LASI', targetIdx);
targetLPSI = getPoints(subj, 'LPSI', targetIdx);
targetRPSI = getPoints(subj, 'RPSI', targetIdx);
targetRASI = rigidBodyFill(baseRPSI, baseLPSI, baseLASI, baseRASI, ...
                           targetRPSI, targetLPSI, targetLASI);
RASI = getPoints(subj, 'RASI');
RASI(targetIdx, :) = targetRASI;
setPoints(subj, 'RASI', RASI);

%% pattern fill (sort off, backward)
subj = 'S03';
[x, y, z, e] = vicon.GetTrajectory(subj, 'LPSI');
baseIdx = find(e, 1, 'last');
baseA = getPoints(subj, 'LPSI', baseIdx);
baseB = getPoints(subj, 'LASI', baseIdx);
v = baseA-baseB;

targetIdx = baseIdx+1:length(e);
targetB = getPoints(subj, 'LASI', targetIdx);
targetA = targetB + v;
LPSI = getPoints(subj, 'LPSI');
LPSI(targetIdx, :) = targetA;
setPoints(subj, 'LPSI', LPSI);

%% pattern fill (sort off, forward)
subj = 'S01';
[x, y, z, e] = vicon.GetTrajectory(subj, 'LPSI');
baseIdx = max(find(e, 1), 2);
baseA = getPoints(subj, 'LPSI', baseIdx);
baseB = getPoints(subj, 'RPSI', baseIdx);
v = baseA-baseB;

targetIdx = 1:baseIdx-1;
targetB = getPoints(subj, 'RPSI', targetIdx);
targetA = targetB + v;
LPSI = getPoints(subj, 'LPSI');
LPSI(targetIdx, :) = targetA;
setPoints(subj, 'LPSI', LPSI);

%% sit
[x, y, z, e] = vicon.GetTrajectory(subj, 'LPSI');
targetIdx = 1:1150;
RASI = getPoints(subj, 'RASI');
LASI = getPoints(subj, 'LASI');
LPSI = getPoints(subj, 'LPSI');

baseIdx = 2662;
baseR = getRotationMatrix(LASI(baseIdx,:), RASI(baseIdx,:), LPSI(baseIdx,:));
v_GCS = (LPSI(baseIdx,:) - LASI(baseIdx,:))'; % 3 x 1 target vector
v_LCS = baseR'*v_GCS; % 3 x 1 vector. components of basis X, Y, Z
                      % v in local coordinate system

for i=targetIdx
    basisX = RASI(i,:)-LASI(i,:);
    basisX = basisX / norm(basisX);
    
    basisZ = baseR(:,3)' - dot(baseR(:,3)', basisX)*basisX;
    basisZ = basisZ / norm(basisZ);
    
    basisY = cross(basisZ, basisX);
    
    estR = [basisX' basisY' basisZ'];
    LPSI(i,:) = (estR*v_LCS)' + LASI(i,:);
end
e(targetIdx) = true;
setPoints(subj, 'LPSI', LPSI, e);