function out = calcCalibAnkleSensorW2PelvisWFromAcc(obj, idx)
	% Calculate yaw offset from accelerometer data
	%
	% Run localization Kalman filter on the pelvis and ankle IMUs
	% Assume the position vector (from origin) of each IMU should point to the
	% same direction, and the yaw angle between vectors is the yaw offset.
	%
	% :param obj: this XsensBody
	% :param idx: index of data to be used
	%
	% :return: out - XsensBody class with adjustment sensor data
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 8/18/18
	
    out = mocapdb.XsensBody();
    
    segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
    
    pos = struct();
    VAR_WIN  = floor(obj.fs*0.25); % NUM_SAMPLES
    ACC_VAR_THRESH = 1;
    A = eye(9, 9);
    A(1:3, 4:6) = eye(3) / obj.fs; 
    A(1:3, 7:9) = - eye(3) / (obj.fs*obj.fs);
    A(4:6, 7:9) = - eye(3) / obj.fs;
    B = zeros(9, 3);
    B(1:3, 1:3) = eye(3) / (obj.fs*obj.fs);
    B(4:6, 1:3) = eye(3) / obj.fs;
    Q = B * 0.5 * eye(3) * B';
    Q(7:9, 7:9) = 0.5;
    H = zeros(6, 9);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);
    R = eye(6) * 1e-4;
    
    for i=1:length(segList)
        n = segList{i};
        
        if sum(size(obj.(n))) ~= 0
            x_post = zeros(9, 1);
            P_post = Q;
            
            freeacc = quatrotate(quatconj(obj.(n).ori(idx, :)), obj.(n).acc(idx, :)) - [0 0 9.81];
            movVarAcc = movingvar(sqrt( sum(freeacc .^2, 2)), VAR_WIN);
            if strcmp(n, 'Pelvis')
                bIsStep = movVarAcc < ACC_VAR_THRESH;
            else
                bIsStep = movVarAcc < ACC_VAR_THRESH;
            end
            
            kfN = length(idx);
            x_prioAll = zeros(kfN, 9); x_postAll = zeros(kfN, 9);
            for kfI = 1:kfN
                x_prio = A * x_post + B * freeacc(kfI, :)';
                P_prio = A * P_post * A' + Q;
            
                if bIsStep(kfI)
                    y = [zeros(3, 1); freeacc(kfI, :)'] - H*x_prio;
                    K = P_prio * H' / (H * P_prio * H' + R);
                    x_post = x_prio + K * y;
                    P_post = (eye(9) - K*H)*P_prio;
                else
                    x_post = x_prio;
                    P_post = P_prio;
                end
                
                x_prioAll(kfI, :) = x_prio';
                x_postAll(kfI, :) = x_post';
            end
            
            pos.(n) = x_postAll;
        end
    end
    
    X = pos.Pelvis(end, 1:2) / norm(pos.Pelvis(end, 1:2));
    for i=1:length(segList)
        n = segList{i};
        
        if sum(size(obj.(n))) ~= 0
            Y = pos.(n)(end, 1:2) / norm(pos.(n)(end, 1:2));
            theta = acosd(dot(X, Y));
            R = [cosd(theta) sind(theta) 0;
                 -sind(theta) cosd(theta) 0;
                 0 0 1];
            out.(n).ori = rotm2quat(R);
        end
    end
end