function out = changePosUnit(obj, newUnit, update)
	% Change position unit of BVH body
	% 
	% Example:
	%      out = obj.changePosUnit('m', false);
	%
	% :param obj: class BVHBody (self)
	% :param newUnit: new unit
	% :param update: If true, update this bvh body, else crease new BVHBody
	%
	% :return: out - BVHBody class whose data only includes the rows in idx
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 11/14/18

    validateattributes(newUnit, {'string', 'char'}, {});
    if nargin <=2, update=false; end
    
    if update
        out = obj;
    else
        out = obj.copy();
    end
    
    if (obj.posUnit == 'mm') & (newUnit == 'm')
        mult = 1.0/1000;
    else
        error('Conversion from %s to %s not yet supported', obj.posUnit, newUnit);
    end
    
    for i=1:length(obj.posList)
        n = obj.posList{i};
        if sum(size(out.(n))) ~= 0
            out.(n) = obj.(n) .* mult;
        end
    end
    out.posUnit = newUnit;
end