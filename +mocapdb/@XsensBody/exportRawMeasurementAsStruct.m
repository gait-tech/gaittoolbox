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
        for j=1:3
            out.(sprintf('%sAcc%c', segAlias{i}, lmap(j))) = t.acc(:,j);
        end
        for j=1:3
            out.(sprintf('%sGyr%c', segAlias{i}, lmap(j))) = t.gyr(:,j);
        end
        for j=1:3
            out.(sprintf('%sMag%c', segAlias{i}, lmap(j))) = t.mag(:,j);
        end
        for j=1:4
            out.(sprintf('%sq%d', segAlias{i}, j)) = t.ori(:,j);
        end
    end
end