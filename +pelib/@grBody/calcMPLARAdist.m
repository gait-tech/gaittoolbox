% ======================================================================
%> @brief Calculate the midpelvis, left ankle, and right ankle distances of grBody
%> @author Luke Sy (UNSW GSBME)
%> @date 15 Mar 2019
%>
%> Example: out = obj.calcMPLARAdist(100)
%> out = struct('MPLA', n x 1, 'MPRA', n x 1, 'LARA', n x 1)
%>
%> @param obj this grBody
%>
%> @retval out struct of distances
% ======================================================================
function out = calcMPLARAdist(obj)    
    out = {};
    out.MPLA = vecnorm(obj.MIDPEL-obj.LTIO, 2, 2);
    out.MPRA = vecnorm(obj.MIDPEL-obj.RTIO, 2, 2);
    out.LARA = vecnorm(obj.LTIO-obj.RTIO, 2, 2);
end