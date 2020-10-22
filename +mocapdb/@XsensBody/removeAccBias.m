function out = removeAccBias(obj, bias)
	% Normalize acceleration
	% 
	% :param obj: this object
    % :type obj: :class:`+mocapdb.@XsensBody`
    % :param bias: XsensBody containing bias information
    % :type obj: :class:`+mocapdb.@XsensBody`
	%
	% :return: out - XsensBody with biased remove
	% :rtype: :class:`+mocapdb.@XsensBody`
    %
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 Jun 10
    
    %% Variable initialization
    out = obj.copy();
    
    %% Calculation
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n)) && ~isempty(obj.(n).acc) ...
           && ~isempty(bias.(n)) && ~isempty(bias.(n).acc)     
            out.(n).acc = obj.(n).acc - bias.(n).acc;
        end
    end
end