function out = adjustFrame(obj, xb1, xb2, orionly)
	% adjust XsensBody frames by xb1.qB * obj.qB * xb2.qB
	% 
	% :param obj: this XsensBody
	% :param qR1: XsensBody or quaternion 1
    % :param qR2: XsensBody or quaternion 2
    % :param orionly: calculate and return only ori values
	%
	% :return: out - XsensBody with adjusted convention.
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 23

    if nargin <= 3
        orionly = false;
    end
    if orionly
        out = mocapdb.XsensBody('srcFileName', obj.srcFileName, 'frame', obj.frame, ...
                                'nSamples', obj.nSamples, 'fs', obj.fs);
    else
        out = obj.copy();
    end
    
    for i=1:length(obj.segList)
        n = obj.segList{i};
        if ~isempty(obj.(n))
            if isa(xb1, 'mocapdb.XsensBody')
                if ~isempty(xb1.(n)), qR1 = xb1.(n).ori;
                else, qR1 = [1 0 0 0]; end
            else
                qR1 = xb1;
            end
            if isa(xb2, 'mocapdb.XsensBody')
                if ~isempty(xb2.(n)), qR2 = xb2.(n).ori;
                else, qR2 = [1 0 0 0]; end
            else
                qR2 = xb2;
            end
            out.(n).ori = quatmultiply(qR1, quatmultiply(obj.(n).ori, qR2));
            
            if ~orionly
                out.(n).acc = quatrotate(qR2, out.(n).acc);
                out.(n).gyr = quatrotate(qR2, out.(n).gyr);
                out.(n).mag = quatrotate(qR2, out.(n).mag);
            end
        end
    end
end