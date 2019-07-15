setup = struct('file', 'S01-Trial-Walk-1', 'algo', "NS2");
ns = extractBetween(setup.algo, 1, 3);
load(sprintf('neura-sparse01/explore-v2/%s-%s-debug.mat', ns, setup.file));
fs = 100;

d = quatmultiply(quatconj(qOri.w__sv.PELV), qOri.w__sv.PELV([2:end end], :));
w = zeros(size(d,1), 3);
for i=2:size(d,1)
    w(i,:) = rot2vec(quat2rotm(d(i,:)));
end

rawAVelPelv = wbodyOri.w__sv.PELV;
oriAVelPelv = w*fs;
pelib.viz.plotXYZ(fs, rawAVelPelv, oriAVelPelv);