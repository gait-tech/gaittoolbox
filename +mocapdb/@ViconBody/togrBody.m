% ======================================================================
%> @brief Export ViconBody class to grBody class
%>
%> @param idx index to be copied to grBody class
%> @param args arguments to be passed to grBody constructor
%>
%> @retval out instance of grBody class.
% ======================================================================
function out = togrBody(obj, idx, args)
    out = pelib.grBody(args{:});
    
    key1 = {'RFEP', 'RFEO', 'RTIO', 'LFEP', 'LFEO', 'LTIO'};
    key2 = {'qRPV', 'qRTH', 'qRSK', 'qLTH', 'qLSK'};
           
    val1 = {'RFEP', 'RFEO', 'RTIO', 'LFEP', 'LFEO', 'LTIO'};
    val2 = {'qRPV', 'qRTH', 'qRSK', 'qLTH', 'qLSK'};
           
    for i=1:length(key1)
        data = (obj.(key1{i}));
        out.(val1{i}) = data(idx,:);
    end
    
    for i=1:length(key2)
        data = (obj.(key2{i}));
        out.(val2{i}) = data(idx,:);
    end
    
    out.MIDPEL = [mean([obj.LFEP(idx,1) obj.RFEP(idx,1)], 2),...
                  mean([obj.LFEP(idx,2) obj.RFEP(idx,2)], 2),...
                  mean([obj.LFEP(idx,3) obj.RFEP(idx,3)], 2)];
              
    out.posUnit = obj.posUnit;
    out.fs = obj.fs;
    out.nSamples = length(idx);
    out.frame = obj.frame;
end