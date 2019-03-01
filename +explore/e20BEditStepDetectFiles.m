% ======================================================================
%> Edit step detect files
% ======================================================================
dir = 'neura-sparse01';
expDir = sprintf('%s/explore-v2', dir);
stepDir = sprintf('%s/step-detect', dir);
ns = "NS2";
algo = "NS2+Aw__sOw__sIw__v+Sav01+M76+C355";
instruction = readtable(sprintf('%s/edit.csv', dir));
dataN = size(instruction, 1);

step = table();
for i = 1:dataN
    n = table2struct(instruction(i, :));
    n.cmd = lower(n.cmd);
    if strcmp(n.cmd, 'load')
        step = readtable(sprintf('%s/%s', stepDir, n.v0));
        name = n.v0(1:strfind(n.v0, '-imuStepDetect')-1);
        load(sprintf('%s/%s-%s-debug.mat', expDir, ns, name));
        load(sprintf('%s/%s-%s-%s.mat', expDir, ns, name, algo));
        idx = allIdx.(cs.initSrc);

        fprintf("Loaded %s\n", n.v0);
    elseif strcmp(n.cmd, 'set') || strcmp(n.cmd, 'clear')
        if strcmp(n.cmd, 'set'), val = 1; 
        else val = 0; end
        if strcmp(lower(n.v0), 'r'), target = 'stepR';
        else target = 'stepL'; end
        if ischar(n.v1), sIdx = str2num(n.v1);
        else sIdx = n.v1; end
        if ischar(n.v2), eIdx = str2num(n.v2);
        else eIdx = n.v2; end
        step.(target)(idx(sIdx):idx(eIdx)) = val;
    elseif strcmp(n.cmd, 'save')
        writetable(step, sprintf('%s/%s', stepDir, n.v0));
        fprintf("Saved %s\n", n.v0);
    end
end