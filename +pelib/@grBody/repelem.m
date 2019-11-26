function out = repelem(obj, n)
% Replicate the elements of grBody
%
% :param obj: class grBody (self)
% :param n: Each element of obj is repeated n (scalar) times. 
%           The new length will be length(obj.element)*n.
%
% :return: out grBody class whose element is repeated n times
% 
% .. Author: - Luke Sy (UNSW GSBME) - 2019/11/26

validateattributes(n, {'numeric'}, {});
    
out = obj.copy();
for i=1:length(out.posList)
    j = out.posList{i};
    out.(j) = repelem(out.(j), n, 1);
end
for i=1:length(out.oriList)
    j = out.oriList{i};
    out.(j) = repelem(out.(j), n, 1);
end
out.nSamples = obj.nSamples*n;

end