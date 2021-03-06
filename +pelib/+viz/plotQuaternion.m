function plotQuaternion(varargin)   
	% Plot s0 s1 s2 s3 of quaternion
	% 
	% :param inputs: array of n x 4's to be plotted
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    sp = [];
    lnames = arrayfun(@inputname, 1:nargin, 'UniformOutput', false);
    
    for i=1:4
        sp(i) = subplot(4,1,i); hold on;
        
        for j=1:nargin
            plot(varargin{j}(:,i));
        end
        
        legend(lnames)
    end
    linkaxes(sp, 'x');
end