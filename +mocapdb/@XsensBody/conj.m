function out = conj(obj)
	% Conjugate of orientation data in XsensBody. Acc, gyr, mag are
	% removed.
	% 
	% :param obj: this XsensBody
	%
	% :return: out - conjugate of XsensBody with adjusted convention.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 23

    out = mocapdb.XsensBody('srcFileName', obj.srcFileName, 'frame', obj.frame, ...
                            'nSamples', obj.nSamples, 'fs', obj.fs);
    
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n))
            out.(n).ori = quatconj(obj.(n).ori);
        end
    end
end