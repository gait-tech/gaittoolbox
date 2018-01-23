% ======================================================================
%> @brief Plot body orientation
%>
%> @param bodys Body instance(s) to be plotted
%> @param parts String(s) of body point(s) to be plotted.
%>
% ======================================================================
function plotOrientation(bodys, parts)
    if ~iscell(bodys)
        bodys = {bodys};
    end
    if ~iscell(parts)
        parts = {parts};
    end
    
    n = length(parts)*3;
    plotIndex = 1;
    
    for i=1:length(parts)
        for j=1:length(bodys)
            data = quat2eul(bodys{j}.(parts{i}));
            if bodys{j}.oriUnit == 'deg'
                data = data * 180 / pi;
            end
            
            t = (1:length(data(:,1)))';

            subplot(n,1,plotIndex); hold on;
            plot(t, data(:,1), ...
                 strcat(bodys{j}.xyzColor{1}, bodys{j}.lnSymbol));
            subplot(n,1,plotIndex+1); hold on;
            plot(t, data(:,2), ...
                 strcat(bodys{j}.xyzColor{2}, bodys{j}.lnSymbol));
            subplot(n,1,plotIndex+2); hold on;
            plot(t, data(:,3), ...
                 strcat(bodys{j}.xyzColor{3}, bodys{j}.lnSymbol));
        end
        
        subplot(n,1,plotIndex);
        title(strcat(parts{i}, ' Z'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        subplot(n,1,plotIndex+1);
        title(strcat(parts{i}, ' Y'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        subplot(n,1,plotIndex+2);
        title(strcat(parts{i}, ' X'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        plotIndex = plotIndex + 3;
    end
    xlabel('Time');
end