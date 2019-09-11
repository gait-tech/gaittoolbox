% ======================================================================
%> @brief generate grBody debug versions
%> @author Luke Sy (UNSW GSBME)
%> @date 27 May 2019
%>
%>
%> @param body grBody
%> @param sensors struct containing sensor data
%> @param debugData debug information
%> @param algo estimator mode
%> @param suffix added to the field name appended to sensorsDebug struct
%>
%> @retval out grBody
% ======================================================================
function [bodyDebug, sensorsDebug] = buildgrBodyDebug(body, sensors, debugData, algo, suffix)
    if nargin <= 4
        suffix = '';
    end    
    bodyDebug = body.copy();
    n2 = body.nSamples*3;
    
    if isstruct(sensors)
        sensorsDebug = struct();
        fn = fieldnames(sensors);
        for i=1:length(fn)
            n = fn{i};
            sensorsDebug.(n) = repelem(sensors.(n), 3, 1);
        end
    else
        sensorsDebug = struct();;
    end
    
    if strcmp(algo, 'vicon')
        bodyDebug = body.copy();
        for i=1:length(bodyDebug.posList)
            n = bodyDebug.posList{i};
            bodyDebug.(n) = repelem(bodyDebug.(n), 3, 1);
        end
        for i=1:length(bodyDebug.oriList)
            n = bodyDebug.oriList{i};
            bodyDebug.(n) = repelem(bodyDebug.(n), 3, 1);
        end
    else
        bodyDebug.MIDPEL = repelem(bodyDebug.MIDPEL, 3, 1);
        bodyDebug.LTIO = repelem(bodyDebug.LTIO, 3, 1);
        bodyDebug.RTIO = repelem(bodyDebug.RTIO, 3, 1);
        bodyDebug.qRPV = repelem(bodyDebug.qRPV, 3, 1);
        bodyDebug.qLSK = repelem(bodyDebug.qLSK, 3, 1);
        bodyDebug.qRSK = repelem(bodyDebug.qRSK, 3, 1);

        if strcmp(algo, 'lieekfv1')
            bodyDebug.MIDPEL(1:3:n2, :) = squeeze(debugData.xhatPri.W_T_PV(1:3,4,:))';
            bodyDebug.MIDPEL(2:3:n2, :) = squeeze(debugData.xhatPos.W_T_PV(1:3,4,:))';
            bodyDebug.LTIO(1:3:n2, :) = squeeze(debugData.xhatPri.W_T_LS(1:3,4,:))';
            bodyDebug.LTIO(2:3:n2, :) = squeeze(debugData.xhatPos.W_T_LS(1:3,4,:))';
            bodyDebug.RTIO(1:3:n2, :) = squeeze(debugData.xhatPri.W_T_RS(1:3,4,:))';
            bodyDebug.RTIO(2:3:n2, :) = squeeze(debugData.xhatPos.W_T_RS(1:3,4,:))';
            bodyDebug.qRPV(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_T_PV(1:3,1:3,:));
            bodyDebug.qRPV(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_T_PV(1:3,1:3,:));
            bodyDebug.qLSK(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_T_LS(1:3,1:3,:));
            bodyDebug.qLSK(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_T_LS(1:3,1:3,:));
            bodyDebug.qRSK(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_T_RS(1:3,1:3,:));
            bodyDebug.qRSK(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_T_RS(1:3,1:3,:));
        elseif strcmp(algo, 'lieekfv2')
            bodyDebug.MIDPEL(1:3:n2, :) = debugData.xhatPri.vec(1:3,:)';
            bodyDebug.MIDPEL(2:3:n2, :) = debugData.xhatPos.vec(1:3,:)';
            bodyDebug.LTIO(1:3:n2, :) = debugData.xhatPri.vec(4:6,:)';
            bodyDebug.LTIO(2:3:n2, :) = debugData.xhatPos.vec(4:6,:)';
            bodyDebug.RTIO(1:3:n2, :) = debugData.xhatPri.vec(7:9,:)';
            bodyDebug.RTIO(2:3:n2, :) = debugData.xhatPos.vec(7:9,:)';
            bodyDebug.qRPV(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_R_PV(1:3,1:3,:));
            bodyDebug.qRPV(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_R_PV(1:3,1:3,:));
            bodyDebug.qLSK(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_R_LS(1:3,1:3,:));
            bodyDebug.qLSK(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_R_LS(1:3,1:3,:));
            bodyDebug.qRSK(1:3:n2, :) = rotm2quat(debugData.xhatPri.W_R_RS(1:3,1:3,:));
            bodyDebug.qRSK(2:3:n2, :) = rotm2quat(debugData.xhatPos.W_R_RS(1:3,1:3,:));
        else
            bodyDebug.MIDPEL(1:3:n2, :) = debugData.predState(:,1:3);
            bodyDebug.MIDPEL(2:3:n2, :) = debugData.zuptState(:,1:3);
            bodyDebug.LTIO(1:3:n2, :) = debugData.predState(:,11:13);
            bodyDebug.LTIO(2:3:n2, :) = debugData.zuptState(:,11:13);
            bodyDebug.RTIO(1:3:n2, :) = debugData.predState(:,21:23);
            bodyDebug.RTIO(2:3:n2, :) = debugData.zuptState(:,21:23);
            bodyDebug.qRPV(1:3:n2, :) = debugData.predState(:,7:10);
            bodyDebug.qRPV(2:3:n2, :) = debugData.zuptState(:,7:10);
            bodyDebug.qLSK(1:3:n2, :) = debugData.predState(:,17:20);
            bodyDebug.qLSK(2:3:n2, :) = debugData.zuptState(:,17:20);
            bodyDebug.qRSK(1:3:n2, :) = debugData.predState(:,27:30);
            bodyDebug.qRSK(2:3:n2, :) = debugData.zuptState(:,27:30);
            
            bodyDebug.qRPV(1, :) = [1 0 0 0];
            bodyDebug.qLSK(1, :) = [1 0 0 0];
            bodyDebug.qRSK(1, :) = [1 0 0 0];
        end

        v = quat2rotm(bodyDebug.qLSK); v = squeeze(v(:,3,:))';
        bodyDebug.LFEO = bodyDebug.LTIO + bodyDebug.getLShankLength(3)*v;
        v = quat2rotm(bodyDebug.qRSK); v = squeeze(v(:,3,:))';
        bodyDebug.RFEO = bodyDebug.RTIO + bodyDebug.getRShankLength(3)*v;
        v = quat2rotm(bodyDebug.qRPV); v = squeeze(v(:,2,:))';
        bodyDebug.LFEP = bodyDebug.MIDPEL + bodyDebug.getPelvisLength(3)/2*v;
        bodyDebug.RFEP = bodyDebug.MIDPEL - bodyDebug.getPelvisLength(3)/2*v;

        v = zeros(3, 3, n2);
        z = (bodyDebug.LFEP-bodyDebug.LFEO)';
        z = z ./ vecnorm(z, 2, 1);
        v(:, 3, :) = reshape(z, 3, 1, []);
        y =  quat2rotm(bodyDebug.qLSK);
        v(:, 2, :) = y(:, 2, :);
        x = cross(v(:, 2, :), v(:, 3, :));
        x = x ./ vecnorm(x, 2, 1);
        v(:, 1, :) =  reshape(x, 3, 1, []);
        bodyDebug.qLTH = rotm2quat(v);

        v = zeros(3, 3, n2);
        z = (bodyDebug.RFEP-bodyDebug.RFEO)';
        z = z ./ vecnorm(z, 2, 1);
        v(:, 3, :) = reshape(z, 3, 1, []);
        y =  quat2rotm(bodyDebug.qRSK);
        v(:, 2, :) = y(:, 2, :);
        x = cross(v(:, 2, :), v(:, 3, :));
        x = x ./ vecnorm(x, 2, 1);
        v(:, 1, :) =  reshape(x, 3, 1, []);
        bodyDebug.qRTH = rotm2quat(v);
        
        if strcmp(algo, 'lieekfv1')
            idx = {{'PELV', 1:3}, {'LANK', 4:6}, {'RANK', 7:9}};
            for i=1:3
                sname = sprintf('%sVel%s', idx{i}{1}, suffix);
                sensorsDebug.(sname)(1:3:n2, :) = debugData.xhatPri.vec(idx{i}{2},:)';
                sensorsDebug.(sname)(2:3:n2, :) = debugData.xhatPos.vec(idx{i}{2},:)';
                sensorsDebug.(sname)(3:3:n2, :) = debugData.xtilde.vec(idx{i}{2},:)';
            end
            
            if size(debugData.xhatPri.vec, 1) > 9
                idx = {{'PELV', 10:12}, {'LANK', 13:15}, {'RANK', 16:18}};
                for i=1:3
                    sname = sprintf('%sAVel%s', idx{i}{1}, suffix);
                    sensorsDebug.(sname)(1:3:n2, :) = debugData.xhatPri.vec(idx{i}{2},:)';
                    sensorsDebug.(sname)(2:3:n2, :) = debugData.xhatPos.vec(idx{i}{2},:)';
                    sensorsDebug.(sname)(3:3:n2, :) = debugData.xtilde.vec(idx{i}{2},:)';
                end
            end
        elseif strcmp(algo, 'lieekfv2')
            idx = 0;
            for i = ["Vel", "AVel"]
                for j = ["PELV", "LANK", "RANK"]
                    k = sprintf("%s%s%s", j, i, suffix);
                    sensorsDebug.(k)(1:3:n2, :) = debugData.xhatPri.vec(idx+(1:3),:)';
                    sensorsDebug.(k)(2:3:n2, :) = debugData.xhatPos.vec(idx+(1:3),:)';
                    sensorsDebug.(k)(3:3:n2, :) = debugData.xtilde.vec(idx+(1:3),:)';
                    idx = idx + 3;
                end
            end
        else
            idx = {{'PELV',  4: 6}, {'LANK', 14:16}, {'RANK', 24:26}};
            for i=1:3
                sname = sprintf('%sVel%s', idx{i}{1}, suffix);
                sensorsDebug.(sname)(1:3:n2, :) = debugData.predState(:,idx{i}{2});
                sensorsDebug.(sname)(2:3:n2, :) = debugData.zuptState(:,idx{i}{2});
                sensorsDebug.(sname)(3:3:n2, :) = debugData.cstrState(:,idx{i}{2});
            end
        end
    end
    bodyDebug.nSamples = bodyDebug.nSamples*3;
    bodyDebug.fs = bodyDebug.fs*3;
end