function [out, bias] = normalizeAcc(obj, acc_norm)
	% Normalize acceleration
	% 
	% :param obj: this object
    % :type obj: :class:`+mocapdb.@XsensBody`
    % :param acc_norm: magnitude of normalized acc. Defaults to 9.81
    % :type acc_norm: Optional, float
	%
	% :return: [out, bias] - XsensBody with normalized acc data and bias data
	% :rtype: [:class:`+mocapdb.@XsensBody`, :class:`+mocapdb.@XsensBody`]
    %
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 Jun 9
    
    if nargin <= 1
        acc_norm = 9.81;
    end
    
    %% Variable initialization
    out = obj.copy();
    bias = obj.copy();
    
    %% Calculation
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if isempty(obj.(n)), continue; end
        
        if ~isempty(obj.(n).acc)
            out.(n).acc = obj.(n).acc ./ vecnorm(obj.(n).acc, 2, 2) .* acc_norm;
            bias.(n).acc = obj.(n).acc - out.(n).acc;
        end
    end
end