% ======================================================================
%> @brief generate c3d file
%> @author Luke Sy (UNSW GSBME)
%> @date 09 Sept 2018
%>
%> Example:
%>      fname = 'test.c3d';
%>      sensors = {'PELVAccX': (n x 1), 'PELVAccY': (n x 1), ... };
%>      refBody = actBody;
%>      
%>      out = obj.exportc3d(fname, sensors, refBody);
%>
%> @param obj grBody (self)
%> @param fname output file name
%> @param sensors [Optional] struct {'label': (n x 1) values to be saved as analog
%>        signals ... } (e.g. raw acc, gyro, magnetometer) 
%> @param refBody [Optional] reference grBody class
%> @param lsteps [Optional] (n x 1) logical where it is true during left foot step detection
%> @param rsteps [Optional] (n x 1) logical where it is true during right foot step detection
%> @param extraMarkers [Optional] extra markers of format struct 
%>                      {'label': (n x 3) position values
%> @param oriMode [Optional] 01: refBody axis on refBody. obj axis on obj.
%>                           02: refBody and obj axis on obj.
%>                           03: refBody and obj axis on refBody.
%> @retval acq handle pointer to new btk c3d
% ======================================================================
function acq = exportc3d(obj, fname, sensors, refBody, lsteps, rsteps, ...
                         extraMarkers, oriMode)
    if nargin <= 2, sensors = struct(); 
    else, validateattributes(sensors, {'struct', 'logical'}, {}); end
    if nargin <= 3, refBody = false; 
    else, validateattributes(refBody, {'pelib.grBody', 'logical'}, {}); end
    
    %% c3d file initializations and metadata
    n = obj.nSamples; fs = obj.fs;
    if nargin <= 4, lsteps = logical(n); end
    if nargin <= 5, rsteps = logical(n); end
    if nargin <= 6, extraMarkers = struct(); end
    if nargin <= 7, oriMode = 1; end
    
    zeroRes = zeros(n, 1);
    acq = btkNewAcquisition(0, n, 0, 1);
    btkSetFrequency(acq, obj.fs);
    btkSetFirstFrame(acq, 1);
    btkSetPointsUnit(acq, 'marker', obj.posUnit);
    btkSetPointsUnit(acq, 'angle', 'deg');
    
    % append meta data
    btkAppendMetaData(acq, 'MANUFACTURER', 'COMPANY', ...
                      btkMetaDataInfo('Char', {'UNSW GSBME'}));
    btkAppendMetaData(acq, 'MANUFACTURER', 'VERSION_LABEL', ...
                      btkMetaDataInfo('Char', {'v1.0-06092018'}));
    desc = struct('MIDPEL', 'Mid pelvis', ...
                  'LFEP', 'Left hip joint center', ... 
                  'LFEO', 'Left knee joint center', ...
                  'LTIO', 'Left ankle joint center', ...
                  'RFEP', 'Right hip joint center', ...
                  'RFEO', 'Right knee joint center', ...
                  'RTIO', 'Right ankle joint center');
              
    %% add markers
    for i = 1:length(obj.posList)
        ptName = obj.posList{i};
        btkAppendPoint(acq, 'marker', ptName, obj.(ptName), zeroRes, ...
                       desc.(ptName));
    end
    if ~islogical(refBody)
        for i = 1:length(refBody.posList)
            ptName = refBody.posList{i};
            btkAppendPoint(acq, 'marker', sprintf('%sRef', refBody.posList{i}), ...
                            refBody.(ptName), zeroRes, desc.(ptName));
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
            if ~islogical(refBody), addAxis(refBody, refBody, acq, 'Ref'); end
        case 2
            addAxis(obj, obj, acq, '');
            if ~islogical(refBody), addAxis(refBody, obj, acq, 'Ref'); end
        case 3
            addAxis(obj, refBody, acq, '');
            if ~islogical(refBody), addAxis(refBody, refBody, acq, 'Ref'); end
    end
    
    % add joint angles
    addAngles(obj, acq, '');
    if ~islogical(refBody)
        addAngles(refBody, acq, 'Ref');
    end
    
    % add analog signals
    analogLabels = fieldnames(sensors);
    for i = 1:length(analogLabels)
        analogName = analogLabels{i};
        btkAppendAnalog(acq, analogName, sensors.(analogName)(1:n));
        
        if analogName(end-3:end-1) == 'Acc'
            btkSetAnalogUnit(acq, analogName, 'm/s^2');
        elseif analogName(end-3:end-1) == 'Gyr'
            btkSetAnalogUnit(acq, analogName, 'rad/s');
        elseif analogName(end-3:end-1) == 'Mag'
            btkSetAnalogUnit(acq, analogName, 'unitless');
        end
    end
    
    % add step detection
    event = lsteps - lsteps([1, 1:end-1]);
    [idx, val] = find(event == 1);
    setEvents(acq, true, idx/fs, 'Left', '', '');
    [idx, val] = find(event == -1);
    setEvents(acq, false, idx/fs, 'Left', '', '');

    event = rsteps - rsteps([1, 1:end-1]);
    [idx, val] = find(event == 1);
    setEvents(acq, true, idx/fs, 'Right', '', '');
    [idx, val] = find(event == -1);
    setEvents(acq, false, idx/fs, 'Right', '', '');
    
    % save c3d file
    btkWriteAcquisition(acq, fname);
end

function setEvents(acq, step, idx, context, subject, desc)
    if step
        label = 'Foot Strike';
        id = 1;
    else
        label = 'Foot Off';
        id = 2;
    end
    for i=1:length(idx)
         btkAppendEvent(acq, label, idx(i), context, subject, desc, id);
    end
end

function addAngles(obj, acq, suffix)
    seq = 'YXZ';
    zeroRes = zeros(obj.nSamples, 1);
    
    [r2 r1 r3] = quat2angle(obj.qRPV, seq); eul = [r1 r2 r3]*180/pi;
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

    [r2 r1 r3] = quat2angle(obj.qLTH, seq); eul = [r1 r2 r3]*180/pi;
    btkAppendPoint(acq, 'angle', sprintf('LThighAngles%s', suffix), eul, ...
                   zeroRes, 'Left thigh angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qRTH, seq); eul = [r1 r2 r3]*180/pi;
    btkAppendPoint(acq, 'angle', sprintf('RThighAngles%s', suffix), eul, ...
                   zeroRes, 'Right thigh angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qLSK, seq); eul = [r1 r2 r3]*180/pi;
    btkAppendPoint(acq, 'angle', sprintf('LShankAngles%s', suffix), eul, ...
                   zeroRes, 'Left shank angles X, Y, Z (seq: YXZ)');
    [r2 r1 r3] = quat2angle(obj.qRSK, seq); eul = [r1 r2 r3]*180/pi;
    btkAppendPoint(acq, 'angle', sprintf('RShankAngles%s', suffix), eul, ...
                   zeroRes, 'Right shank angles X, Y, Z (seq: YXZ)');
end

function addAxis(body1, body2, acq, suffix)
    d = norm(body2.MIDPEL(1,:) - body2.LFEP(1,:))*0.5;
    
    pair = struct('qRPV', 'MIDPEL', 'qLTH', 'LFEO', 'qRTH', 'RFEO', ...
                  'qLSK', 'LTIO', 'qRSK', 'RTIO');
              
    labels = fieldnames(pair);
    for i=1:length(labels)
        v = labels{i}; k = pair.(v);
        R = quat2rotm(body1.(v));
        btkAppendPoint(acq, 'marker', sprintf('%sX%s', v, suffix), ...
                       body2.(k)+d*squeeze(R(:,1,:))' );
        btkAppendPoint(acq, 'marker', sprintf('%sY%s', v, suffix), ...
                       body2.(k)+d*squeeze(R(:,2,:))' );
        btkAppendPoint(acq, 'marker', sprintf('%sZ%s', v, suffix), ...
                       body2.(k)+d*squeeze(R(:,3,:))' );
    end
end