%% Recalculate neura error from saved states

dir = 'neura-sparse01/explore-v2';
dataList = ls(sprintf('%s/neura-*+C*.mat', dir));

expression = 'neura-(?<subj>\w+)-(?<act>[-a-zA-z0-9]+)-NS1\+A(?<acc>\w+)O(?<ori>\w+)I(?<initSrc>\w+)\+Sav01\+M(?<meas>\d+)\+C(?<cstr>\d+)\.mat*';

clear results;
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
    csActBodyRel = csActBody.changeRefFrame('MIDPEL');
    
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
    csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    results0 = estBody2.diffRMSE(csActBody2);
        
    results0.name = name;
    results0.label = sprintf("NS1+A%sO%sI%s+Sav01+M%02d+C%03d", ...
        n.acc, n.ori, n.initSrc, str2num(n.meas), str2num(n.cstr));
    results0.runtime = runtime;
    results(i) = results0;
    display(sprintf("%s/%s-%s", dir, name, results0.label));
end

results = struct2table(results);
save(sprintf("%s/results.mat", dir), 'results')

function label = getLabel(setup)
    if setup.accData == 'v'
        if setup.accDataNoise == 0 
            aD = 'v';
        else
            aD = strrep(sprintf('v%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData;
    end
    label = sprintf('N%s%s%s+M%02d+C%03d', aD, setup.oriData, setup.initSrc, ...
        setup.applyMeas, setup.applyCstr);
end