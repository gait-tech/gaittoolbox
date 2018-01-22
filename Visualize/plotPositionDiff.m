% ======================================================================
%> @brief plot body position difference (body2 - body1)
%>
%> @param body1 Body instance basis
%> @param body2 Body instance(s) to be compared from the basis
%> @param parts String(s) of body point(s) to be plotted.
%> @param symbol2 plot symbol to be used for the graph
%>
% ======================================================================
function plotPositionDiff(body1, body2, parts, symbol2)
    n = (nargin-2)*3;
    plotIndex = 1;
    
    for i=3:nargin
        pts = body.(vargin{i});
        t = length(pts(:,1));
        subplot(n,1,plotIndex);
        plot(t, pts(:,1), strcat('r', symbol));
        subplot(n,1,plotIndex+2);
        plot(t, pts(:,2), strcat('g', symbol));
        subplot(n,1,plotIndex+3);
        plot(t, pts(:,3), strcat('b', symbol));
        
        plotIndex = plotIndex + 3;
    end
end