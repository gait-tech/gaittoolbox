% motion list
mot = struct('file', 'S01-Trial-Walk-1', ...
             'algo', "NS2+lieekfv1+Aw__vOw__vIw__v+Sav03+P003+M000+C000");

dataSfname = sprintf('neura-sparse01/imu/%s', mot.file);
ns = extractBetween(mot.algo, 1, 3);

load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, mot.file));
load(sprintf('neura-sparse01/explore-v2/%s-%s-%s.mat', ns, mot.file, mot.algo));

if cs.initSrc == 'w__v'
    aLabel = 'w__v';
    vb = W__viconBody;
elseif cs.initSrc == 'v__v'
    aLabel = 'v__v';
    vb = V__viconBody;
else
    aLabel = 'w__x';
    vb = W__xsensBody;
end
if ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
   strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
    eLabel = strcat(cs.accData, cs.initSrc(end));
else
    eLabel = cs.accData;
end
idx = allIdx.(cs.initSrc);
    
avel = wbodyOri.(aLabel);
csQOri = qOri.(aLabel);

% avel.LTIB2 = quatrotate(quatconj(csQOri.LTIB), avel.LTIB);
% v = vb.LFEO - vb.LTIO;
% 
% vel.LFEO2 = vel.LTIO + cross(avel.LTIB2, v(:, :));
% clf; pelib.viz.plotXYZ(100, vel.LFEO, vel.LFEO2)

vHmkbasisT = {}; vecIdx = 1;
W_avel = {}; avelIdx = 1;

for body = {vb, estBody}
    body1 = body{1};
    vel = body1.calcJointVel({'MIDPEL', 'LFEO', 'RFEO', 'LFEP', 'RFEP', 'LTIO', 'RTIO'});
    W_avel{avelIdx} = body1.calcSegAngVel({'qRPV', 'qLSK', 'qRSK'}, 'W');
    avelIdx = avelIdx + 1;
    
    W_R_LT = quat2rotm(body1.qLTH);
    W_R_RT = quat2rotm(body1.qRTH);
    vLHmLK = vel.LFEP - vel.LFEO;
    vRHmRK = vel.RFEP - vel.RFEO;
    
    vHmkbasisT{vecIdx} = [dot(vLHmLK, squeeze(W_R_LT(:,1,:))', 2), ...
                          dot(vLHmLK, squeeze(W_R_LT(:,2,:))', 2), ...
                          dot(vLHmLK, squeeze(W_R_LT(:,3,:))', 2)];
    vecIdx = vecIdx + 1;
    
    vHmkbasisT{vecIdx} = [dot(vRHmRK, squeeze(W_R_RT(:,1,:))', 2), ...
                          dot(vRHmRK, squeeze(W_R_RT(:,2,:))', 2), ...
                          dot(vRHmRK, squeeze(W_R_RT(:,3,:))', 2)];
    vecIdx = vecIdx + 1;
end

vLHmLKbasisLTAct = vHmkbasisT{1};
vRHmRKbasisRTAct = vHmkbasisT{2};
vLHmLKbasisLTEst = vHmkbasisT{3};
vRHmRKbasisRTEst = vHmkbasisT{4};

updateFigureContents('vHmK');
pelib.viz.plotXYZ(100, vLHmLKbasisLTAct, vRHmRKbasisRTAct, ...
                       vLHmLKbasisLTEst, vRHmRKbasisRTEst);

WavelActPV = W_avel{1}.qRPV;
WavelEstPV = W_avel{2}.qRPV;
WavelActLS = W_avel{1}.qLSK;
WavelEstLS = W_avel{2}.qLSK;
WavelActRS = W_avel{1}.qRSK;
WavelEstRS = W_avel{2}.qRSK;

updateFigureContents('avel check');
t=1:length(idx);
pelib.viz.plotXYZ(100, WavelActPV, WavelEstPV, WavelActLS, WavelEstLS, ...
                  WavelActRS, WavelEstRS);