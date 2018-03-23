% ======================================================================
%> @brief Compare the state across time
%> @author Luke Sy
%>
%> 
%> @param estStateData estimated state data (n_samples x n_state)
%> @param actStateData actual (basis) state data (n_samples x n_state)
%> @param stateIdx compare state number stateIdx: stateIdx+2
% ======================================================================

function plotStateComparison(estStateData, actStateData, stateIdx)
    plotIndex = 1;
    target = stateIdx:stateIdx+2;
    sp = [];
    n = length(estStateData.predState(:,1))*3;
    
    for i=target
        sp(plotIndex) = subplot(length(target), 1, plotIndex); hold on;
%         t = [1:3:n 2:3:n 3:3:n];
%         y = [estStateData.predState(:,i); ...
%              estStateData.zuptState(:,i); ...
%              estStateData..constrainedState(:,i)];
        scatter(1:3:n, estStateData.predState(:,i), '.r');
        scatter(2:3:n, estStateData.zuptStateL(:,i), '<g');
        scatter(2:3:n, estStateData.zuptStateR(:,i), '>g');
        scatter(3:3:n, estStateData.cstrState(:,i), '.b');
        
        actStateData2 = repelem(actStateData(:,i), 3);
        scatter(1:n, actStateData2, '.k');
        
        legend('pred', 'zuptL', 'zuptR', 'cpknee', 'actual');
        plotIndex = plotIndex + 1;
    end
    linkaxes(sp, 'x');

end