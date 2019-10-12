function out = exportRawMeasurementAsStruct(obj, seg, segAlias)
	% Export raw xsens measurements as struct
	%
	% Example:
	%      seg = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
	%
	%      segAlias = {'PELV', 'LANK', 'RANK'};
	%
	%      out = obj.exportRawMeasurementAsStruct(seg, segAlias);
	%
	% :param obj: class XsensBody (self)
	% :param seg: cell array of sensors to be exported 
	% :param segAlias: cell array of sensor aliases
	% :return: out - {'PELVAccX': n x 1, PELVAccY: n x 1, ... }
	%          for Acc, Gyr, Mag on the X, Y, Z axes
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/9/18


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