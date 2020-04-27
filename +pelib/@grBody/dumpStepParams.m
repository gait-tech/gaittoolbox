function dumpStepParams(obj1, obj2, intervals, fname)
	% Dump each individual stride length and gait speed to a file
	%
	% :param obj1: grBody 1 (self)
	% :param obj2: grBody 2 (other, basis)
	% :param intervals: interval when total distance is calculated
    % :param fname: output step data to files (fname-<posFieldname>.csv)
    %
	%
	% .. Author: - Luke Sy (UNSW GSBME) 2020/04/14

    seq = 'YXZ';
    
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
        T = table();
        % time
        time = (intervals2.(i)(2:end)-intervals2.(i)(1:end-1)+1)/obj1.fs;
        if ~isempty(time)
            % totaltime = (intervals2.(i)(end)-intervals2.(i)(1)+1)/obj1.fs;
            T.("StartTime_s") = intervals2.(i)(1:end-1)/obj1.fs;
            T.("EndTime_s") = intervals2.(i)(2:end)/obj1.fs;
            T.("Duration_s") = time;

            %% obj1 TTD
            % get position at intervals
            pos1 = obj1.(i)(intervals2.(i), 1:2);
            % get distance between intervals
            dist1 = vecnorm(pos1(2:end,:)-pos1(1:end-1,:), 2, 2);
            T.("StrideLengthEst_m") = dist1;
            T.("GaitSpeedEst_mpersec") = dist1./time;

            %% act TTD
            % get position at intervals
            pos2 = obj2.(i)(intervals2.(i), 1:2);
            % get distance between intervals
            dist2 = vecnorm(pos2(2:end,:)-pos2(1:end-1,:), 2, 2);
            % step length
            T.("StrideLengthAct_m") = dist2;
            % gait speed SL
            T.("GaitSpeedAct_mpersec") = dist2./time;
        end
        writetable(T, sprintf('%s-%s.csv',fname,i));
    end
end