function plotXYZ(fs, varargin)   
	% Plot x y z of inputs
	% 
	% :param fs: sampling frequency
	% :param inputs: array of n x 3's to be plotted
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    sp = [];
    lnames = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
    
    for i=1:3
        sp(i) = subplot(3,1,i); hold on;
        
        for j=1:nargin-1
            n = length(varargin{j}(:,i));
            t = (1:n)/fs;
            plot(t, varargin{j}(:,i));
        end
        
        legend(lnames)
    end
    linkaxes(sp, 'x');
end