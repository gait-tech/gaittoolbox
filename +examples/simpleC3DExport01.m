%% Sample to generate c3d export
% This scripts generate the c3d export from the output of runSample01.m
% 
% Note: run this script only after executing runSample01.m
dir = 'data/sample';
expDir = sprintf('%s/output', dir);
ns = "NS2";
n = struct('subj', 'S02', 'act', 'Trial-Walk-1', ...
           'startFrame', 1120, 'endFrame', 2074);
       
name = sprintf("%s-%s-%s", ns, n.subj, n.act);
label = "NS2+ckf+Aw__sOw__sIw__v+Sav03+P001+M076+C355";

load(sprintf("%s/%s-debug.mat", expDir, name));

sIdx = max(allIdx.w__v(1), allIdx.w__x(1));
eIdx = min(allIdx.w__v(end), allIdx.w__x(end));
nSamples = eIdx - sIdx + 1;
viconIdx0 = find(allIdx.w__v==sIdx,1):find(allIdx.w__v==eIdx,1);
xsensIdx0 = find(allIdx.w__x==sIdx,1):find(allIdx.w__x==eIdx,1);

csActBody = W__viconBody.getSubset(viconIdx0);
xsensBody = W__xsensBody.getSubset(xsensIdx0);

targetname = sprintf('%s/%s-viconvsxsens.c3d', expDir, name);
xsensBody.exportc3d(targetname, struct(), csActBody);

load(sprintf('%s/%s-%s.mat', expDir, name, label));
targetname = sprintf('%s/%s-%s.c3d', expDir, name, label);
estBody.exportc3d(targetname, struct(), csActBody);
