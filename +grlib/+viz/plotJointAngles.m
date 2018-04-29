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
    
    for i=1:length(parts)
        for j=1:nBody
            switch lower(parts{i})
                case 'lhip'
                    data = bodys{j}.calcJointAnglesLHip();
                case 'rhip'
                    data = bodys{j}.calcJointAnglesRHip();
                case 'lknee'
                    data = bodys{j}.calcJointAnglesLKnee();
                case 'rknee'
                    data = bodys{j}.calcJointAnglesRKnee();                    
            end
            
            data = data * 180 / pi;
            
            t = (1:length(data(:,1)))';
            t = t / bodys{j}.fs;
            
            subplot(nParts, 1, plotIndex); hold on;
            plot(t, data(:,1), ...
                 strcat(bodys{j}.xyzColor{1}, bodys{j}.lnSymbol));
            subplot(nParts, 1, plotIndex+1); hold on;
            plot(t, data(:,2), ...
                 strcat(bodys{j}.xyzColor{2}, bodys{j}.lnSymbol));
            subplot(nParts, 1, plotIndex+2); hold on;
            plot(t, data(:,3), ...
                 strcat(bodys{j}.xyzColor{3}, bodys{j}.lnSymbol));
        end
        
        if nBody == 2
            
        end
        
        for j=1:3
            subplot(nParts, 1, plotIndex+j-1);
            if nBody == 2
                title(sprintf('%s %c (rmse=%.6f %c)', parts{i}, ...
                              axisName{j}, 0, bodys{1}.posUnit));
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