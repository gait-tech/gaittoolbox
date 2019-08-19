% ======================================================================
%> @brief export raw xsens measurements as struct
%> @author Luke Sy (UNSW GSBME)
%> @date 09 Sept 2018
%>
%> Example:
%>      seg = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
%>      segAlias = {'PELV', 'LANK', 'RANK'};
%>      out = obj.exportRawMeasurementAsStruct(seg, segAlias);
%>
%>      out returns {'PELVAccX': n x 1, PELVAccY: n x 1, ... }
%>          for Acc, Gyr, Mag on the X, Y, Z axes
%>
%> @param obj class XsensBody (self)
%> @param seg cell array of sensors to be exported 
%> @param segAlias cell array of sensor aliases
% ======================================================================
function out = exportRawMeasurementAsStruct(obj, seg, segAlias)
    validateattributes(seg, {'cell'}, {});
    validateattributes(segAlias, {'cell'}, {});
    if length(seg) ~= length(segAlias)
        error('seg and segAlias have different lengths');
    end
    
    g = [0 0 9.81];
    
    out = struct();
    lmap = 'XYZ';
    for i=1:length(seg)
        t = obj.(seg{i});
        out.(sprintf('B_%sAcc', segAlias{i})) = t.acc;
        out.(sprintf('B_%sGyr', segAlias{i})) = t.gyr;
        out.(sprintf('B_%sMag', segAlias{i})) = t.mag;
        out.(sprintf('W_%sAcc', segAlias{i})) = quatrotate(quatconj(t.ori), t.acc);
        out.(sprintf('W_%sGyr', segAlias{i})) = quatrotate(quatconj(t.ori), t.gyr);
        out.(sprintf('W_%sMag', segAlias{i})) = quatrotate(quatconj(t.ori), t.mag);
        
        for j=1:4
            out.(sprintf('%sq%d', segAlias{i}, j)) = t.ori(:,j);
        end
    end
end