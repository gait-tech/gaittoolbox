function out = calcTTDandStepParams(obj1, obj2, intervals, baseStruct)
	% Returns the total travelled distance, stride length, and gait speed
	%
	% :param obj1: grBody 1 (self)
	% :param obj2: grBody 2 (other, basis)
	% :param intervals: interval when total distance is calculated
    % :param baseStruct: [Optional] append input struct if defined
    %
	% :return: out - struct with the error report
	%
	% .. Author: - Luke Sy (UNSW GSBME) 2020/04/09

    seq = 'YXZ';
    if nargin <= 3, out = struct;
    else, out = baseStruct; end
    
    posFields = ["MIDPEL", "LTIO", "RTIO"];
    intervals2 = struct();
    for i=posFields
        buf = intervals.(i);
        buf([1, end]) = true;
        event = buf - buf([1, 1:end-1]);
        [idx, val] = find(event == 1);
        intervals2.(i) = idx;
    end
    
    for i=posFields
        % time
        time = (intervals2.(i)(2:end)-intervals2.(i)(1:end-1)+1)/obj1.fs;
        totaltime = (intervals2.(i)(end)-intervals2.(i)(1)+1)/obj1.fs;
        
        %% obj1 TTD
        % get position at intervals
        pos1 = obj1.(i)(intervals2.(i), 1:2);
        % get distance between intervals
        dist1 = vecnorm(pos1(2:end,:)-pos1(1:end-1,:), 2, 2);
        % add distance
        ttd1 = sum(dist1);
        % step length
        out.(sprintf("%sStrideLengthEstMean",i)) = mean(dist1);
        out.(sprintf("%sStrideLengthEstStd",i)) = std(dist1);
        % gait speed SL
        out.(sprintf("%sGaitSpeedStrLenEstMean",i)) = mean(dist1./time);
        out.(sprintf("%sGaitSpeedStrLenEstStd",i)) = std(dist1./time);
        
        %% act TTD
        % get position at intervals
        pos2 = obj2.(i)(intervals2.(i), 1:2);
        % get distance between intervals
        dist2 = vecnorm(pos2(2:end,:)-pos2(1:end-1,:), 2, 2);
        % add distance
        ttd2 = sum(dist2);
        % step length
        out.(sprintf("%sStrideLengthActMean",i)) = mean(dist2);
        out.(sprintf("%sStrideLengthActStd",i)) = std(dist2);
        % gait speed SL
        out.(sprintf("%sGaitSpeedStrLenActMean",i)) = mean(dist2./time);
        out.(sprintf("%sGaitSpeedStrLenActStd",i)) = std(dist2./time);
        
        % ratio
        out.(sprintf("%sTTDEst",i)) = ttd1;
        out.(sprintf("%sTTDAct",i)) = ttd2;
        out.(sprintf("%sTTDPE",i)) = abs(ttd1/ttd2-1.0);
        % gait speed TTD
        out.(sprintf("%sGaitSpeedTTDEst",i)) = ttd1/totaltime;
        out.(sprintf("%sGaitSpeedTTDAct",i)) = ttd2/totaltime;
    end
end