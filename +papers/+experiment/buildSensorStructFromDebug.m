% ======================================================================
%> @brief generate sensor struct debug versions
%> @author Luke Sy (UNSW GSBME)
%> @date 27 May 2019
%>
%>
%> @param sensors struct containing sensor data
%> @param state debug information
%> @param algo estimator mode
%>
%> @retval out grBody
% ======================================================================
function sensors = buildSensorStructFromDebug(sensors, state, state2, algo, suffix)
    if nargin <= 4
        suffix = '';
    end
    
    if strcmp(algo, 'lieekfv1')
        idx = {{'PELV', 1:3}, {'LANK', 4:6}, {'RANK', 7:9}};
        for i=1:3
            sname = sprintf('%sVel%s', idx{i}{1}, suffix);
            sensors.(sname) = state.vec(idx{i}{2}, :)';
        end
        
        if size(state.vec, 1) > 9
            idx = {{'PELV', 10:12}, {'LANK', 13:15}, {'RANK', 16:18}};
            for i=1:3
                sname = sprintf('%sAVel%s', idx{i}{1}, suffix);
                sensors.(sname) = state.vec(idx{i}{2}, :)';
            end
        end
        
        idx = {{'PV',  1: 3,  4: 6, 19:21, 28:30}, ...
               {'LS',  7: 9, 10:12, 22:24, 31:33}, ...
               {'RS', 13:15, 16:18, 25:27, 34:36}};
        if size(state.vec, 1) > 9
            idx2 = {'Pos', 'Ori', 'Vel', 'AVel'};
        else
            idx2 = {'Pos', 'Ori', 'Vel'};
        end
        [nState, ~, nSamples] = size(state2.Ptilde);
        for i=1:3
            buf = zeros(nSamples, nState);
            for j=1:nSamples
                buf(j,:) = diag(state2.Ptilde(:,:,j));
            end
            for j=1:length(idx2)
            	sname = sprintf('P%s%s%s', idx2{j}, idx{i}{1}, suffix);
            	sensors.(sname) = buf(:,idx{i}{j+1});
                
                sname2 = sprintf('vPos%s%s%s', idx2{j}, idx{i}{1}, suffix);
                sensors.(sname2) = state2.measUptPos(idx{i}{j+1}, :)';
                
                sname3 = sprintf('vTilde%s%s%s', idx2{j}, idx{i}{1}, suffix);
                sensors.(sname3) = state2.measUptTilde(idx{i}{j+1}, :)';
            end
        end
        
        idx = {{'PV',  1: 3}, {'LS',  4: 6}, {'RS', 7: 9}};
        for i=1:3
            sname = sprintf('TAcc%s%s', idx{i}{1}, suffix);
            sensors.(sname) = state2.u(idx{i}{2}, :)';
        end
    elseif strcmp(algo, 'lieekfv2')
        idx = 9;
        for i = ["Vel", "AVel"]
            for j = ["PELV", "LANK", "RANK"]
                k = sprintf("%s%s%s", j, i, suffix);
                sensors.(k) = state.vec(idx+(1:3), :)';
                idx = idx + 3;
            end
        end

        idx = 0;         
        [nState, ~, nSamples] = size(state2.Ptilde);
        buf = zeros(nSamples, nState);
        for i=1:nSamples
            buf(i,:) = diag(state2.Ptilde(:,:,i));
        end
        for i = ["Ori", "Pos", "Vel", "AVel"]
            for j = ["PV", "LS", "RS"]
                k = sprintf("P%s%s%s", i, j, suffix);
                sensors.(k) = buf(:,idx+(1:3));
                
                k = sprintf("vPos%s%s%s", i, j, suffix);
                sensors.(k) = state2.measUptPos(idx+(1:3), :)';
                
                k = sprintf("vTilde%s%s%s", i, j, suffix);
                sensors.(k) = state2.measUptTilde(idx+(1:3), :)';
                idx = idx + 3;
            end
        end
        
        idx = 0;
        for i = ["PV", "LS", "RS"]
            k = sprintf("uAcc%s%s", i, suffix);
            sensors.(k) = state2.u(idx+(1:3), :)';
            idx = idx + 3;
        end
    else
        idx = {{'PELV',  4: 6}, {'LANK', 14:16}, {'RANK', 24:26}};
        for i=1:3
            sname = sprintf('%sVel%s', idx{i}{1}, suffix);
            sensors.(sname) = state(:, idx{i}{2});
        end
    end
end