function out = togrBody(obj, idx, args)
	% Export BVHBody class to grBody class
	%
	% :param idx: index to be copied to grBody class
	% :param args: arguments to be passed to grBody constructor
	%
	% :return: out - instance of grBody class.
	%
	% .. Author: - Luke Sy (UNSW GSBME)
	
    out = pelib.grBody(args{:});
    
    key1 = {'RightUpLeg', 'RightLeg', 'RightFoot', 'RightToe', ...
            'LeftUpLeg', 'LeftLeg', 'LeftFoot', 'LeftToe'};
    key2 = {'qHips', 'qRightUpLeg', 'qRightLeg', 'qRightFoot', ...
            'qLeftUpLeg', 'qLeftLeg', 'qLeftFoot'};
           
    val1 = {'RFEP', 'RFEO', 'RTIO', 'RTOE', 'LFEP', 'LFEO', 'LTIO', 'LTOE'};
    val2 = {'qRPV', 'qRTH', 'qRSK', 'qRFT', 'qLTH', 'qLSK', 'qLFT'};
           
    for i=1:length(key1)
        data = (obj.(key1{i}));
        out.(val1{i}) = data(idx,:);
    end
    
    % convert TCD ori convention to biomechanical convention
%     qTCD2BM = rotm2quat([0 -1 0; -1 0 0; 0 0 -1]);
    for i=1:length(key2)
        data = (obj.(key2{i}));
%         out.(val2{i}) = quatmultiply(data(idx,:), qTCD2BM);
        out.(val2{i}) = data(idx,:);
    end
    
    out.MIDPEL = [mean([obj.LeftUpLeg(idx,1) obj.RightUpLeg(idx,1)], 2),...
                  mean([obj.LeftUpLeg(idx,2) obj.RightUpLeg(idx,2)], 2),...
                  mean([obj.LeftUpLeg(idx,3) obj.RightUpLeg(idx,3)], 2)];
              
    out.posUnit = obj.posUnit;
    out.nSamples = length(idx);
    out.frame = obj.frame;
end