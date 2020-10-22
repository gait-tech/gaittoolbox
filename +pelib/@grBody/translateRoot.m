function out = translateRoot(obj, pos)
	% Translate body
	%
	% :param obj: grBody (self)
	% :param pos: 3d translation
    %
	% :return: out - new grBody in world frame
	% %rtype: :class:`+pelib.@grBody`
    %
	% .. Author: - Luke Sy (UNSW GSBME) 2020 Jun 10
    
    
    if ~strcmp(obj.frame, 'world')
        error('The source frame %s is not supported', obj.frame);    
    end
    
    out = obj.copy();
    
    posList = obj.posList;
    for i=1:length(posList)
        if(~isempty(obj.(posList{i})))
            out.(posList{i}) = obj.(posList{i}) + pos;
        end
    end
end