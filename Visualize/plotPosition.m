% ======================================================================
%> @brief plot body position
%>
%> @param body Body instance to be plotted
%> @param parts String(s) of body point(s) to be plotted.
%> @param symbol Plot symbol to be used for the graph
%>
% ======================================================================
function plotPosition(body, parts, symbol)
    if ~iscell(parts)
        parts = {parts};
    end
    if nargin <= 2
        symbol = '-';
    end
    
    n = length(parts)*3;
    plotIndex = 1;
    
    for i=3:nargin
        data = body.(parts{i});
        t = length(data(:,1));
        
        subplot(n,1,plotIndex);
        plot(t, data(:,1), strcat('r', symbol));
        subplot(n,1,plotIndex+2);
        plot(t, data(:,2), strcat('g', symbol));
        subplot(n,1,plotIndex+3);
        plot(t, data(:,3), strcat('b', symbol));
        
        plotIndex = plotIndex + 3;
    end
end