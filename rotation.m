classdef rotation < handle
    
    % INTERFACE ROTATION CLASS:
    %
    % This class represents an interface class for the rotation hyerarchy
    %
    % Author: Gabriele Ligorio
    % Version: 1.0
    % Date: 04-11-2015

    
    properties (Access = protected)
        value;
        len;
    end
    
    methods (Abstract)
        
        % Virtual methods
        rslt = getValue(obj);
        rslt = setValue(obj,in,pos);
        rslt = getN(obj);
        rslt = inv(obj);
        rslt = unit(obj); 
        rslt = rep(obj,N);
        rslt = average(obj,win);
        
    end
    
end

