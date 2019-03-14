%% Recalculate neura error from saved states
dir0 = 'neura-sparse01';
dir = sprintf('%s/explore-v2', dir0);
ns = 'NS2';
clear results;
rIdx = 1;
% results = table2struct(results);
% rIdx = 757;

%% file list vicon vs xsens comparison
dataList = readtable(sprintf('%s/data-list-v2.csv', dir0));
dataN = size(dataList, 1);

for i = 1:dataN
    n = table2struct(dataList(i, :));
    name = sprintf("%s-%s-%s", ns, n.subj, n.act);
    load(sprintf("%s/%s-debug.mat", dir, name));
    nSamples = min(W__viconBody.nSamples, W__xsensBody.nSamples);
    csActBody = W__viconBody.getSubset(1:nSamples);
    estBody = W__xsensBody.getSubset(1:nSamples);
    
    
    csActBodyRel = csActBody.changeRefFrame('MIDPEL');
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
    csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    results0 = estBody2.diffRMSEandMean(csActBody2);
        
    results0.name = name;
    results0.label = sprintf("%s+viconvsxsens", ns);
    results0.runtime = 0;
    results(rIdx) = results0;
    rIdx = rIdx + 1;
    fprintf("%s/%s-%s\n", dir, name, results0.label);
end

%% folder .mat rerun
dataList = ls(sprintf('%s/%s-*+C*.mat', dir, ns));
expression = 'NS.-(?<subj>\w+)-(?<act>[-a-zA-z0-9]+)-NS.\+A(?<acc>\w+)O(?<ori>\w+)I(?<initSrc>\w+)\+Sav01\+M(?<meas>\d+)\+C(?<cstr>\d+)\.mat*';

for i = 1:size(dataList, 1)
    load(sprintf('%s/%s', dir, dataList(i, :)));
    n = regexp(dataList(i, :), expression, 'names');
    if strcmp(string(n.subj), ""), continue; end
    
    name = sprintf('%s-%s-%s', 'neura', n.subj, n.act);
    load(sprintf("%s/%s-debug.mat", dir, name));
    
    if n.initSrc == 'w__v'
        csActBody = W__viconBody;
    elseif n.initSrc == 'v__v'
        csActBody = V__viconBody;
    else
        csActBody = W__xsensBody;
    end
    nSamples = min(csActBody.nSamples, estBody.nSamples);
    csActBody = csActBody.getSubset(1:nSamples);
    estBody = estBody.getSubset(1:nSamples);
    
    csActBodyRel = csActBody.changeRefFrame('MIDPEL');
    
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
    csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    results0 = estBody2.diffRMSEandMean(csActBody2);
        
    results0.name = name;
    results0.label = sprintf("%s+A%sO%sI%s+Sav01+M%02d+C%03d", ...
        ns, n.acc, n.ori, n.initSrc, str2num(n.meas), str2num(n.cstr));
    results0.runtime = runtime;
    results(rIdx) = results0;
    rIdx = rIdx + 1;
    display(sprintf("%s/%s-%s", dir, name, results0.label));
end

results = struct2table(results);
save(sprintf("%s/results.mat", dir), 'results')

function label = getLabel(ns, setup)
    if setup.accData == 'v'
        if setup.accDataNoise == 0 
            aD = 'v';
        else
            aD = strrep(sprintf('v%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData;
    end
    label = sprintf("%s+A%sO%sI%s+S%s+M%02d+C%03d", ns, aD, setup.oriData, setup.initSrc, ...
        setup.stepDetection, setup.applyMeas, setup.applyCstr);
end