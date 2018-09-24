markerData = struct('Names', containers.Map, 'Pos', []);
idx = 3;
for i=1:length(viconBody.posList)
    n = viconBody.posList{i};
    markerData.Names(n) = idx:idx+2;
    markerData.Pos(:,idx:idx+2) = viconBody.(n)*1000;
    idx = idx + 3;
end

markerSet = {'LFEP', 'MIDPEL'; 'RFEP', 'MIDPEL';
             'LFEO', 'LFEP'; 'RFEO', 'RFEP';
             'LFEO', 'LTIO'; 'RFEO', 'RTIO'};
markerSetColour = {'g', 'g', 'r', 'b', 'r', 'b'};



pelib.viz.animViconMarkersV2('markerData', markerData, ...
    'markerSet', markerSet, 'markerSetColour', markerSetColour);