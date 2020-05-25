function out = calcTTD(obj1, obj2, intervals, baseStruct)
	% Returns the RMSE and Mean difference between grBody1 and grBody2
	%
	% :param obj1: grBody 1 (self)
	% :param obj2: grBody 2 (other, basis)
	% :param intervals: interval when total distance is calculated
    % :param baseStruct: [Optional] append input struct if defined
    %
	% :return: out - struct with the difference of pos and ori parameters
	%
	% .. Author: - Luke Sy (UNSW GSBME)

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
        %% obj1 TTD
        % get position at intervals
        pos1 = obj1.(i)(intervals2.(i), 1:2);
        % get distance between intervals
        dist1 = vecnorm(pos1(2:end,:)-pos1(1:end-1,:), 2, 2);
        % add distance
        ttd1 = sum(dist1);
        
        %% act TTD
        % get position at intervals
        pos2 = obj2.(i)(intervals2.(i), 1:2);
        % get distance between intervals
        dist2 = vecnorm(pos2(2:end,:)-pos2(1:end-1,:), 2, 2);
        % add distance
        ttd2 = sum(dist2);
        
        % ratio
        out.(sprintf("%sTTDEst",i)) = ttd1;
        out.(sprintf("%sTTDAct",i)) = ttd2;
        out.(sprintf("%sTTDPE",i)) = abs(ttd1/ttd2-1.0);
    end
end