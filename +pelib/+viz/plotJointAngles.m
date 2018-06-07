% ======================================================================
%> @brief Plot joint angles
%>
%> @param bodys Body instance(s) to be plotted
%> @param parts String(s) of body point(s) to be plotted (LHip, RHip, LKnee, RKnee).
%>
% ======================================================================
function plotJointAngles(bodys, parts)
    if ~iscell(bodys)
        bodys = {bodys};
    end
    if ~iscell(parts)
        parts = {parts};
    end
    
    nParts = length(parts)*3;
    nBody = length(bodys);
    plotIndex = 1;
    axisName = {'x', 'y', 'z'};
    
    data = {};
    for i=1:length(parts)
        for j=1:nBody
            switch lower(parts{i})
                case 'lhip'
                    data{j} = bodys{j}.calcJointAnglesLHip();
                case 'rhip'
                    data{j} = bodys{j}.calcJointAnglesRHip();
                case 'lknee'
                    data{j} = bodys{j}.calcJointAnglesLKnee();
                case 'rknee'
                    data{j} = bodys{j}.calcJointAnglesRKnee();                    
            end
            
            if bodys{1}.oriUnit == 'deg'
                dataBuf = data{j} * 180 / pi;
            else
                dataBuf = data{j};
            end
            
            t = (1:length(dataBuf(:,1)))';
            t = t / bodys{j}.fs;
            
            subplot(nParts, 1, plotIndex); hold on;
            plot(t, dataBuf(:,1), ...
                 strcat(bodys{j}.xyzColor{1}, bodys{j}.lnSymbol));
            subplot(nParts, 1, plotIndex+1); hold on;
            plot(t, dataBuf(:,2), ...
                 strcat(bodys{j}.xyzColor{2}, bodys{j}.lnSymbol));
            subplot(nParts, 1, plotIndex+2); hold on;
            plot(t, dataBuf(:,3), ...
                 strcat(bodys{j}.xyzColor{3}, bodys{j}.lnSymbol));
        end
        
        if nBody == 2
            errVal = grlib.rmse(grlib.anglediff(data{1}, data{2}));
        end
        
        for j=1:3
            subplot(nParts, 1, plotIndex+j-1);
            if nBody == 2
                title(sprintf('%s %c (rmse=%.6f %s)', parts{i}, ...
                              axisName{j}, errVal(j), bodys{1}.oriUnit));
            else
                title(sprintf('%s %c', parts{i}, axisName{j}));
            end
            legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
            ylabel('Degrees');
        end
        
        plotIndex = plotIndex + 3;
    end
    xlabel('Time');
end