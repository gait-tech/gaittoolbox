% ======================================================================
%> @brief Plot dots based on mean MU and covariance SIGMA
%> @author Luke Sy
%>
%> Plot the lower body
%> 
%> @param MU mean
%> @param SIGMA convariance matrix
%> @param N number of dots to represent the covariance cloud
%> @param setup plot setup
%> @retval p plot object
% ======================================================================
function p = plotMeanCovSamples(MU, SIGMA, N, setup)    
    if nargin <= 2, N = 500; end
    if nargin <= 3, setup = 'b.'; end
    hold on;
    
    try
        SIGMA2 = (SIGMA + SIGMA.') / 2;
        samples = mvnrnd(MU, SIGMA2, N);
        scatter3(samples(:, 1), samples(:, 2), samples(:, 3), setup);
    end
end