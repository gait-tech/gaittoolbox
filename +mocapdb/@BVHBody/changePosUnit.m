% ======================================================================
%> @brief change position unit of BVH body
%> @author Luke Sy (UNSW GSBME)
%> @date 14 Nov 2018
%>
%> Example:
%>      out = obj.changePosUnit('m', false);
%>
%> @param obj class BVHBody (self)
%> @param newUnit new unit
%> @param update If true, update this bvh body, else crease new BVHBody
%>
%> @retval out BVHBody class whose data only includes the rows in idx
% ======================================================================
function out = changePosUnit(obj, newUnit, update)
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