function out = togrBody(obj, idx, args)
	% Export ViconBody class to grBody class
	%
	% :param idx: index to be copied to grBody class
	% :param args: arguments to be passed to grBody constructor
	%
	% :return: out - instance of grBody class.
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    out = pelib.grBody(args{:});
    
    key1 = {'RFEP', 'RFEO', 'RTIO', 'RTOE', 'LFEP', 'LFEO', 'LTIO', 'LTOE'};
    key2 = {'qRPV', 'qRTH', 'qRSK', 'qLTH', 'qLSK', 'qLFT', 'qRFT'};
           
    val1 = {'RFEP', 'RFEO', 'RTIO', 'RTOE', 'LFEP', 'LFEO', 'LTIO', 'LTOE'};
    val2 = {'qRPV', 'qRTH', 'qRSK', 'qLTH', 'qLSK', 'qLFT', 'qRFT'};
           
    for i=1:length(key1)
        data = (obj.(key1{i}));
        if(~isempty(data))
            out.(val1{i}) = data(idx,:);
        end
    end
    
    for i=1:length(key2)
        data = (obj.(key2{i}));
        if(~isempty(data))
            out.(val2{i}) = data(idx,:);
        end
    end
    
    out.MIDPEL = [mean([obj.LFEP(idx,1) obj.RFEP(idx,1)], 2),...
                  mean([obj.LFEP(idx,2) obj.RFEP(idx,2)], 2),...
                  mean([obj.LFEP(idx,3) obj.RFEP(idx,3)], 2)];
              
    out.posUnit = obj.posUnit;
    out.fs = obj.fs;
    out.nSamples = length(idx);
    out.frame = obj.frame;
    out.ftStartIndex = obj.ftStartIndex;
    out.ftEndIndex = obj.ftEndIndex;
end