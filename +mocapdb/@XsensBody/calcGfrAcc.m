function out = calcGfrAcc(obj)
	% Calculate ground frame acceleration (w/out gravity) 
	%
	% :return: out - struct with gfrAcc of each body segment
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 23 May 2020
    
    out = struct();
    for i=obj.segList
        if ~isempty(obj.(i{1}))
            out.(i{1}) = quatrotate(quatconj(obj.(i{1}).ori), obj.(i{1}).acc) ...
                            - [0 0 9.81];
        end
    end
end