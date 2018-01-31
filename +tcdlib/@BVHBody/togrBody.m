% ======================================================================
%> @brief Export BVHBody class to grBody class
%>
%> @param idx index to be copied to grBody class
%> @param args arguments to be passed to grBody constructor
%>
%> @retval out instance of grBody class.
% ======================================================================
function out = togrBody(obj, idx, args)
    out = grBody(args{:});
    
    key1 = {'Hips', 'RightUpLeg', 'RightLeg', 'RightFoot', ...
            'LeftUpLeg', 'LeftLeg', 'LeftFoot'};
    key2 = {'qHips', 'qRightUpLeg', 'qRightLeg', 'qRightFoot', ...
            'qLeftUpLeg', 'qLeftLeg', 'qLeftFoot'};
           
    val1 = {'MIDPEL', 'RFEP', 'RFEO', 'RTIO', 'LFEP', 'LFEO', 'LTIO'};
    val2 = {'qRPV', 'qRTH', 'qRSK', 'qRFT', 'qLTH', 'qLSK', 'qLFT'};
           
    for i=1:length(key1)
        data = (obj.(key1{i}));
        out.(val1{i}) = data(idx,:);
    end
    for i=1:length(key2)
        data = (obj.(key2{i}));
        out.(val2{i}) = data(idx,:);
    end
    out.posUnit = obj.posUnit;
    out.nSamples = length(idx);
end