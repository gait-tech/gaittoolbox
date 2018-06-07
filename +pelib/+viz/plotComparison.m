% ======================================================================
%> @brief Plot comparison between data1 and data2
%> @author Luke Sy
%>
%> @param data1 array of n x m to be plotted
%> @param data2 array of n x m to be plotted
%> @param fs sampling frequency
% ======================================================================

function plotComparison(data1, data2, fs)
    if nargin <= 2
        fs = 1;
    end
    
    sp = [];
    lnames = arrayfun(@inputname, 1:2, 'UniformOutput', false);
    
    [n, m] = size(data1);
    t = (1:n)/fs;
    for i=1:m
        sp(i) = subplot(m,1,i); hold on;
        
        plot(t, data1(:,i), t, data2(:,i));
        
        lnames2 = cellfun(@(x) sprintf('%s-%d', x,i), lnames, ...
                           'UniformOutput', false);
        legend(lnames2)
    end
    linkaxes(sp, 'x');
end