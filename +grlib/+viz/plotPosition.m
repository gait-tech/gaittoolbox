% ======================================================================
%> @brief Plot body position
%>
%> @param bodys Body instance(s) to be plotted
%> @param parts String(s) of body point(s) to be plotted.
%>
% ======================================================================
function plotPosition(bodys, parts)
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
            data = bodys{j}.(parts{i});
            t = (1:length(data(:,1)))';
            t = t / bodys{j}.fs;
            
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
        title(strcat(parts{i}, ' x'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.posUnit);
        
        subplot(n,1,plotIndex+1);
        title(strcat(parts{i}, ' y'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.posUnit);
        
        subplot(n,1,plotIndex+2);
        title(strcat(parts{i}, ' z'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.posUnit);
        
        plotIndex = plotIndex + 3;
    end
    xlabel('Time');
end