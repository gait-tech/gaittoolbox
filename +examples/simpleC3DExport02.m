%% Sample to generate c3d export
% This scripts generate the c3d export from ViconBody and XsensBody
% 

dir = 'data/sample';
n = struct('subj', 'S02', 'act', 'Trial-Walk-1');

dataV = mocapdb.ViconBody.loadViconMat( ...
            sprintf('%s/vicon/%s-%s.mat', dir, n.subj, n.act));           
dataX = mocapdb.BVHBody.loadXsensBVHFile( ...
            sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act), 'mm');
fs = dataV.fs;

V__viconBody = dataV.togrBody(1:dataV.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
V__viconBody.exportc3d(sprintf('%s-%s-vicon.c3d',  n.subj, n.act));

qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
dataX = dataX.toWorldFrame(qXsensV2W);
W__xsensBody = dataX.togrBody(1:dataX.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
W__xsensBody.exportc3d(sprintf('%s-%s-xsens.c3d',  n.subj, n.act));
