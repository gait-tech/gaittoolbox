function acq = exportc3d(obj, fname, sensors, refBody, lsteps, rsteps, ...
                         extraMarkers, oriMode, spevents)
	% Generate c3d file
	%
	% Example:
	%      fname = 'test.c3d';
	%      sensors = {'PELVAccX': (n x 1), 'PELVAccY': (n x 1), ... };
	%      refBody = actBody;
	%      
	%      out = obj.exportc3d(fname, sensors, refBody);
	%
	% :param obj: grBody (self)
    % :type obj: :class:`+pelib.@grBody`
	% :param fname: output file name
	% :param sensors: {'label': (n x 1) values to be saved as analog
	%        signals ... } (e.g. raw acc, gyro, magnetometer) 
    % :type sensors: Optional, struct.
	% :param refBody: reference grBody class
    % :type refBody: Optional, :class:`+pelib.@grBody` or array of :class:`+pelib.@grBody`
	% :param lsteps: (n x 1) logical where it is true during left foot step detection
    % :type lsteps: Optional, logical array
	% :param rsteps: (n x 1) logical where it is true during right foot step detection
    % :type rsteps: Optional, logical array
	% :param extraMarkers: [Optional] extra markers of format struct 
	%                      {'label': (n x 3) position values
	% :param oriMode: 01 - refBody axis on refBody. obj axis on obj.
	%                 02 - refBody and obj axis on obj.
	%                 03 - refBody and obj axis on refBody.
    % :param oriMode: Optional, integer. Default to 1.
	% :param spevents: (n x 1) logical when general event happens
    % :param spevents: Optional, logical array
	% :return: acq - handle pointer to new btk c3d
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 9/09/18
	
    if nargin <= 2, sensors = struct(); 
    else, validateattributes(sensors, {'struct', 'logical'}, {}); end
    if nargin <= 3, refBody = false; 
    else, validateattributes(refBody, {'pelib.grBody', 'logical', 'cell'}, {}); end
    if islogical(refBody)
        refBody = {};
    elseif iscell(refBody)
    else
        refBody = {refBody};
    end
    
    %% c3d file initializations and metadata
    n = obj.nSamples; fs = obj.fs;
    if nargin <= 4, lsteps = false(n, 1); end
    if nargin <= 5, rsteps = false(n, 1); end
    if nargin <= 6, extraMarkers = struct(); end
    if nargin <= 7, oriMode = 1; end
    if nargin <= 8, spevents = false(n, 1); end
    
    zeroRes = zeros(n, 1);
    acq = btkNewAcquisition(0, n, 0, 1);
    btkSetFrequency(acq, obj.fs);
    btkSetFirstFrame(acq, 1);
    btkSetPointsUnit(acq, 'marker', obj.posUnit);
    btkSetPointsUnit(acq, 'angle', 'deg');
    btkSetPointsUnit(acq, 'scalar', 'n/a (see description)');
    
    % append meta data
    btkAppendMetaData(acq, 'MANUFACTURER', 'COMPANY', ...
                      btkMetaDataInfo('Char', {'UNSW GSBME'}));
    btkAppendMetaData(acq, 'MANUFACTURER', 'VERSION_LABEL', ...
                      btkMetaDataInfo('Char', {'v1.0-06092018'}));
    desc = struct('MIDPEL', 'Mid pelvis', ...
                  'LFEP', 'Left hip joint center', ... 
                  'LFEO', 'Left knee joint center', ...
                  'LTIO', 'Left ankle joint center', ...
                  'LTOE', 'Left toe marker', ...
                  'RFEP', 'Right hip joint center', ...
                  'RFEO', 'Right knee joint center', ...
                  'RTIO', 'Right ankle joint center', ...
                  'RTOE', 'Left toe marker' );
              
    %% add markers
    for i = 1:length(obj.posList)
        ptName = obj.posList{i};
        if(~isempty(obj.(ptName)))
            btkAppendPoint(acq, 'marker', ptName, obj.(ptName), zeroRes, ...
                           desc.(ptName));
        end
        % only put in reference point if the main point exists
        for j = 1:length(refBody)
            if(~isempty(refBody{j}.(ptName)))
                btkAppendPoint(acq, 'marker', sprintf('%s%s', ptName, getPostfix(j)), ...
                               refBody{j}.(ptName), zeroRes, desc.(ptName));
            end
        end
    end

    eMarkerLabels = fieldnames(extraMarkers);
    for i = 1:length(eMarkerLabels)
        ptName = eMarkerLabels{i};
        btkAppendPoint(acq, 'marker', ptName, extraMarkers.(ptName));
    end
    
    % add axis
    switch (oriMode)
        case 1
            addAxis(obj, obj, acq, '');
            for j=1:length(refBody)
                addAxis(refBody{j}, refBody{j}, acq, getPostfix(j)); 
            end
        case 2
            addAxis(obj, obj, acq, '');
            for j=1:length(refBody)
                addAxis(refBody{j}, obj, acq, getPostfix(j)); 
            end
        case 3
            addAxis(obj, refBody, acq, '');
            for j=1:length(refBody)
                addAxis(refBody{j}, refBody{j}, acq, getPostfix(j)); 
            end
    end
    
    % add joint angles
    addAngles(obj, acq, '');
    for j=1:length(refBody)
        addAngles(refBody{j}, acq, getPostfix(j));
    end
    
    % add sensor signals
    analogLabels = sort(fieldnames(sensors));
    for i = 1:length(analogLabels)
        analogName = analogLabels{i};
        if size(sensors.(analogName), 2) == 1
            btkAppendAnalog(acq, analogName, sensors.(analogName)(1:n));
            
            if contains(analogName(end-3:end), 'Acc')
                btkSetAnalogUnit(acq, analogName, 'm/s^2');
            elseif contains(analogName(end-3:end), 'Gyr')
                btkSetAnalogUnit(acq, analogName, 'rad/s');
            elseif contains(analogName(end-3:end), 'Mag')
                btkSetAnalogUnit(acq, analogName, 'unitless');
            end
        elseif size(sensors.(analogName), 2) == 3
            if contains(analogName(end-3:end), 'Acc')
                desc = 'm/s^2';
            elseif contains(analogName(end-3:end), 'Gyr')
                desc = 'rad/s';
            elseif contains(analogName(end-3:end), 'Mag')
                desc = 'unitless';
            else
                desc = '';
            end
            btkAppendPoint(acq, 'scalar', analogName, sensors.(analogName), ...
                    zeroRes, desc);
        else
            warning('Cannot append raw measurement %s', analogName); 
        end

    end
    
    % add step detection
    event = lsteps - lsteps([1, 1:end-1]);
    [idx, val] = find(event == 1);
    setEvents(acq, true, (idx-1)/fs, 'Left', '', '', 'Foot Strike');
    [idx, val] = find(event == -1);
    setEvents(acq, false, (idx-1)/fs, 'Left', '', '', 'Foot Off');

    event = rsteps - rsteps([1, 1:end-1]);
    [idx, val] = find(event == 1);
    setEvents(acq, true, (idx-1)/fs, 'Right', '', '', 'Foot Strike');
    [idx, val] = find(event == -1);
    setEvents(acq, false, (idx-1)/fs, 'Right', '', '', 'Foot Off');
    
    event = spevents- spevents([1, 1:end-1]);
    [idx, val] = find(event == 1);
    setEvents(acq, true, (idx-1)/fs, 'General', '', '', 'Start');
    [idx, val] = find(event == -1);
    setEvents(acq, false, (idx-1)/fs, 'General', '', '', 'End');
    
    % save c3d file
    btkWriteAcquisition(acq, char(fname));
end

function setEvents(acq, step, idx, context, subject, desc, label)
    if step
        id = 1;
    else
        id = 2;
    end
    for i=1:length(idx)
         btkAppendEvent(acq, label, idx(i), context, subject, desc, id);
    end
end

function addAngles(obj, acq, suffix)
    seq = 'YXZ';
    zeroRes = zeros(obj.nSamples, 1);
    
    [r2 r1 r3] = quat2angle(obj.qRPV, seq); eul = rad2deg([r1 r2 r3]);
    btkAppendPoint(acq, 'angle', sprintf('PelvisAngles%s', suffix), eul, ...
                    zeroRes, 'Pelvis angles X, Y, Z');
    btkAppendPoint(acq, 'angle', sprintf('LHipAngles%s', suffix), ...
                   obj.calcJointAnglesLHip()*180/pi, ...
                   zeroRes, 'Left hip angles X, Y, Z');
    btkAppendPoint(acq, 'angle', sprintf('RHipAngles%s', suffix), ...
                   obj.calcJointAnglesRHip()*180/pi, ...
                   zeroRes, 'Right hip angles X, Y, Z');
    btkAppendPoint(acq, 'angle', sprintf('LKneeAngles%s', suffix), ...
                   obj.calcJointAnglesLKnee()*180/pi, ...
                   zeroRes, 'Left knee angles X, Y, Z');
    btkAppendPoint(acq, 'angle', sprintf('RKneeAngles%s', suffix), ...
                   obj.calcJointAnglesRKnee()*180/pi, ...
                   zeroRes, 'Right knee angles X, Y, Z');
    if ~isempty(obj.qLFT)
        btkAppendPoint(acq, 'angle', sprintf('LAnkleAngles%s', suffix), ...
                   obj.calcJointAnglesLAnkle()*180/pi, ...
                   zeroRes, 'Left ankle angles X, Y, Z');
    end
    if ~isempty(obj.qRFT)
        btkAppendPoint(acq, 'angle', sprintf('RAnkleAngles%s', suffix), ...
                   obj.calcJointAnglesRAnkle()*180/pi, ...
                   zeroRes, 'Right ankle angles X, Y, Z');
    end

    [r2 r1 r3] = quat2angle(obj.qLTH, seq); eul = rad2deg([r1 r2 r3]);
    btkAppendPoint(acq, 'angle', sprintf('LThighAngles%s', suffix), eul, ...
                   zeroRes, 'Left thigh angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qRTH, seq); eul = rad2deg([r1 r2 r3]);
    btkAppendPoint(acq, 'angle', sprintf('RThighAngles%s', suffix), eul, ...
                   zeroRes, 'Right thigh angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qLSK, seq); eul = rad2deg([r1 r2 r3]);
    btkAppendPoint(acq, 'angle', sprintf('LShankAngles%s', suffix), eul, ...
                   zeroRes, 'Left shank angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qRSK, seq); eul = rad2deg([r1 r2 r3]);
    btkAppendPoint(acq, 'angle', sprintf('RShankAngles%s', suffix), eul, ...
                   zeroRes, 'Right shank angles X, Y, Z (seq: YXZ)');
end

function addAxis(body1, body2, acq, suffix)
    d = norm(body2.MIDPEL(1,:) - body2.LFEP(1,:))*0.5;
    
    pair = struct('qRPV', 'MIDPEL', 'qLTH', 'LFEO', 'qRTH', 'RFEO', ...
                  'qLSK', 'LTIO', 'qRSK', 'RTIO', ...
                  'qLFT', 'LTOE', 'qRFT', 'RTOE');
              
    labels = fieldnames(pair);
    for i=1:length(labels)
        v = labels{i}; k = pair.(v);
        if(~isempty(body1.(v)) && ~isempty(body2.(k)))
            R = quat2rotm(body1.(v));
            btkAppendPoint(acq, 'marker', sprintf('%sX%s', v, suffix), ...
                           body2.(k)+d*squeeze(R(:,1,:))' );
            btkAppendPoint(acq, 'marker', sprintf('%sY%s', v, suffix), ...
                           body2.(k)+d*squeeze(R(:,2,:))' );
            btkAppendPoint(acq, 'marker', sprintf('%sZ%s', v, suffix), ...
                           body2.(k)+d*squeeze(R(:,3,:))' );
        end
    end
end

function out = getPostfix(idx)
    if idx == 1
        out = 'Ref';
    else
        out = sprintf('%d', idx);
    end
end