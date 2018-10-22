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
    
    out = struct();
    lmap = 'XYZ';
    for i=1:length(seg)
        t = obj.(seg{i});
        out.(sprintf('%sAcc', segAlias{i})) = t.acc;
        out.(sprintf('%sGyr', segAlias{i})) = t.gyr;
        out.(sprintf('%sMag', segAlias{i})) = t.mag;
        out.(sprintf('%sAccWorldFrame', segAlias{i})) = quatrotate(quatconj(t.ori), t.acc);
        out.(sprintf('%sGyrWorldFrame', segAlias{i})) = quatrotate(quatconj(t.ori), t.gyr);
        out.(sprintf('%sMagWorldFrame', segAlias{i})) = quatrotate(quatconj(t.ori), t.mag);
        
        for j=1:4
            out.(sprintf('%sq%d', segAlias{i}, j)) = t.ori(:,j);
        end
    end
end