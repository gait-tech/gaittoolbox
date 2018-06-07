% ======================================================================
%> @brief Plot x y z of inputs
%> @author Luke Sy
%>
%> @param inputs array of n x 3's to be plotted
% ======================================================================

function plotXYZ(varargin)   
    sp = [];
    lnames = arrayfun(@inputname, 1:nargin, 'UniformOutput', false);
    
    for i=1:3
        sp(i) = subplot(3,1,i); hold on;
        
        for j=1:nargin
            plot(varargin{j}(:,i));
        end
        
        legend(lnames)
    end
    linkaxes(sp, 'x');
end