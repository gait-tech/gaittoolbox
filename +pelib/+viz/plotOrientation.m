function plotOrientation(bodys, parts)
	% Plot body orientation
	%
	% :param bodys: Body instance(s) to be plotted
	% :param parts: String(s) of body point(s) to be plotted.
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    if ~iscell(bodys)
        bodys = {bodys};
    end
    if ~iscell(parts)
        parts = {parts};
    end
    
    n = length(parts)*3;
    plotIndex = 1;
    axs = []; axIdx = 1;
    
    for i=1:length(parts)
        for j=1:length(bodys)
            data = quat2eul(bodys{j}.(parts{i}));
            if bodys{j}.oriUnit == 'deg'
                data = data * 180 / pi;
            end
            
            t = (1:length(data(:,1)))';

            axs(axIdx) = subplot(n,1,plotIndex); hold on; axIdx = axIdx + 1;
            plot(t, data(:,1), ...
                 strcat(bodys{j}.xyzColor{1}, bodys{j}.lnSymbol));
            axs(axIdx) = subplot(n,1,plotIndex+1); hold on; axIdx = axIdx + 1;
            plot(t, data(:,2), ...
                 strcat(bodys{j}.xyzColor{2}, bodys{j}.lnSymbol));
            axs(axIdx) = subplot(n,1,plotIndex+2); hold on; axIdx = axIdx + 1;
            plot(t, data(:,3), ...
                 strcat(bodys{j}.xyzColor{3}, bodys{j}.lnSymbol));
        end
        
        axs(axIdx) = subplot(n,1,plotIndex); axIdx = axIdx + 1;
        title(strcat(parts{i}, ' Z'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        axs(axIdx) = subplot(n,1,plotIndex+1); axIdx = axIdx + 1;
        title(strcat(parts{i}, ' Y'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        axs(axIdx) = subplot(n,1,plotIndex+2); axIdx = axIdx + 1;
        title(strcat(parts{i}, ' X'));
        legend(cellfun(@(x) x.name, bodys, 'UniformOutput', false));
        ylabel(bodys{1}.oriUnit);
        
        plotIndex = plotIndex + 3;
    end
    xlabel('Time');
    linkaxes(axs, 'x');
end