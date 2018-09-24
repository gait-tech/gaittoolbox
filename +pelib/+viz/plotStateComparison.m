% ======================================================================
%> @brief Compare the state across time
%> @author Luke Sy
%>
%> 
%> @param estStateData estimated state data (n_samples x n_state)
%> @param actStateData actual (basis) state data (n_samples x n_state)
%> @param stateIdx compare state number stateIdx: stateIdx+2
% ======================================================================

function plotStateComparison(estStateData, actStateData, stateIdx, fs)
    plotIndex = 1;
    target = stateIdx:stateIdx+2;
    sp = [];
    n = length(estStateData.predState(:,1))*3;
    fsAdj = fs*3;
    
    for i=target
        sp(plotIndex) = subplot(length(target), 1, plotIndex); hold on;
%         t = [1:3:n 2:3:n 3:3:n];
%         y = [estStateData.predState(:,i); ...
%              estStateData.zuptState(:,i); ...
%              estStateData..constrainedState(:,i)];
        scatter((1:3:n)/fsAdj, estStateData.predState(:,i), '.r');
        zuptState = estStateData.zuptState(:,i);
        zuptStateL = zuptState;
        zuptStateL(~estStateData.zuptStateL(:,1)) = nan;
        zuptStateR = zuptState;
        zuptStateR(~estStateData.zuptStateR(:,1)) = nan;
        scatter((2:3:n)/fsAdj, zuptState(:,1), '.g');
        scatter((2:3:n)/fsAdj, zuptStateL(:,1), '<g');
        scatter((2:3:n)/fsAdj, zuptStateR(:,1), '>g');
        cstrState = estStateData.cstrState(:,i);
%         cstrStateU = sum(estStateData.cstrStateU, 2) > 0;
%         cstrState(~cstrStateU) = nan;
        scatter((3:3:n)/fsAdj, cstrState, '.b');
        
        actStateData2 = repelem(actStateData(:,i), 3);
        scatter((1:n)/fsAdj, actStateData2, '.k');
        
        legend('pred', 'zupt', 'zuptL', 'zuptR', 'cstr', 'actual');
        plotIndex = plotIndex + 1;
    end
    linkaxes(sp, 'x');

end