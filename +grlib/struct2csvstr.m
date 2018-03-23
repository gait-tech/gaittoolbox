function out = struct2csvstr(s, showcolname)
    if nargin <= 1
        showcolname = false;
    end
    
    data = struct2table(s);
    outBuf = {};
    outIdx = 1;
    
    if showcolname
        outBuf{outIdx} = strjoin(data.Properties.VariableNames, ',,,');
        outIdx = outIdx + 1;
    end
    
    n = size(data.Variables); n = n(1);
    for i=1:n
        outBuf{outIdx} = strjoin(fillmissing(string(data{i,:}), 'constant', "nan"), ',');
        outIdx = outIdx + 1;
    end
    
    out = strjoin(string(outBuf), '\n');
end