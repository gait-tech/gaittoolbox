fs=100;
optionTIB = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40BA5', 'R_LowLeg', '00B40C35', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

% optionANK = struct('Pelvis', '00B40B91', ...
%     'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
%     'L_LowLeg', '00B40C49', 'R_LowLeg', '00B40C4A', ...
%     'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
optionANK = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');
xyz = {'X', 'Y', 'Z'};
dir = 'neura/imu/';
% sensorlist = {horzcat(arrayfun(@(x) sprintf('S01-Trial-%03d', x), 0:21, 'UniformOutput', false), ...
%     arrayfun(@(x) sprintf('S02-Trial-%03d', x), 0:20, 'UniformOutput', false), ...
%     arrayfun(@(x) sprintf('S03-Trial-%03d', x), 0:21, 'UniformOutput', false) )};
sensorlist = {horzcat(arrayfun(@(x) sprintf('S03-Trial-%03d', x), 0:21, 'UniformOutput', false) )};
sensorlist = sensorlist{1};

for i=1:length(sensorlist)
    fnameS = sprintf('%s%s', dir, sensorlist{i});
    
    sensTIB = mocapdb.XsensBody.loadMTExport(fnameS, optionTIB);
    sensANK = mocapdb.XsensBody.loadMTExport(fnameS, optionANK);

    d1 = quatrotate(quatconj(sensANK.L_LowLeg.ori), sensANK.L_LowLeg.mag);
    d2 = quatrotate(quatconj(sensTIB.L_LowLeg.ori), sensTIB.L_LowLeg.mag);
    d3 = quatrotate(quatconj(sensANK.R_LowLeg.ori), sensANK.R_LowLeg.mag);
    d4 = quatrotate(quatconj(sensTIB.R_LowLeg.ori), sensTIB.R_LowLeg.mag);

    for j=1:3
        subplot(3,1,j)   
        plot([d1(:,j), d2(:,j), d3(:,j), d4(:,j)]);
        legend({'left ankle', 'left tibia', 'right ankle', 'right tibia'});
        title(sprintf('%s %c', fnameS, xyz{j}));
    end
    
    saveas(gcf, sprintf('explore_output/%s.png', sensorlist{i}));
end