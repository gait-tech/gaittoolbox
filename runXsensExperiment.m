% ======================================================================
%> @brief Run experiment on an instance of Xsens data
%> @author Luke Sy
%> 
%>
%> Setup parameters:
%> - label: data instance name (e.g. s1-acting1 or s2-walking1)
%> - est: filter type to be used.
%>      - ekfv3: grlib.est.kf_3_kmus_v3
%> - accData: acceleration data to be used
%>      - raw: raw xsens measurements
%>      - sim: simulated xsens measurements
%> - oriData: orientation data to be used
%>      - x: xsens
%> - applyMeas: measurement configuration number
%> - applyCstr: constraint configuration number
%> - sigmaQAcc: Q acceleration sigma (variance)
%> - P: initial P matrix
%>
%> @param fnameBVH bvh filename or loaded tcdlib.BVHBody 
%>              (e.g. totalcapture/vicon/s1/acting1_BlenderZXY_YmZ.bvh)
%> @param fnameRAW raw xsens measurement filename or loaded tcdlib.XsensBody
%>              (e.g. totalcapture/gyroMag/s1/Acting1_Xsens_AuxFields.sensors)
%> @param name name of the experiment
%> @param setups list of experiment parameters (struct) to be run. see
%> details above
%> @param savedir filepath to save .mat output/debug files (optional)
% ======================================================================
function results = runXsensExperiment(fnameBVH, fnameRaw, name, setups, savedir)
    %% Inputs and Input Check
    validateattributes(fnameBVH, {'char', 'tcdlib.BVHBody'}, {});
    validateattributes(fnameRaw, {'char', 'tcdlib.XsensBody'}, {});
    if nargin <= 4
        savedir = ''
    end
    
    %% Initialization   
    % Load video orientation and position for each body segment
    if ischar(fnameBVH)
        dataBVH = tcdlib.BVHBody.loadBVHFile(fnameBVH, 'mm');
    else
        dataBVH = fnameBVH;
    end     
    
    % Load sensor orientation and position for each body segment
    if ischar(fnameRaw) 
        dataRaw = tcdlib.XsensBody.loadSensorFile(fnameRaw);
    else
        dataRaw = fnameRaw;
    end
    
    % Initialize other variables
    fs = 100;
    qV2W = rotm2quat([1 0 0; 0 0 -1; 0 1 0]);
    
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'raw', 'accDataNoise', 0.0, 'oriData', 'x', ...
        'stepDetectWindow', 0.25, 'stepDetectThreshold', 1, ...
        'meas', 0, 'cstr', 0, ...
        'sigmaQAcc', 0.5, 'P', 100);
    setupN = length(setups);
    setupDefaultFN = fieldnames(setupDefault);
        
    for i=1:setupN
        for j=1:length(setupDefaultFN)
            if ~isfield(setups{i}, setupDefaultFN{j})
                setups{i}.(setupDefaultFN{j}) = setupDefault.(setupDefaultFN{j});
            end
        end
    end
    
    % convert vicon to world frame
    nSamples = min(dataBVH.nSamples, dataRaw.nSamples);
    dataBVH = dataBVH.toWorldFrame(qV2W);
    key = {'Pelvis', 'L_UpLeg', 'R_UpLeg',...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    val = {'Hips', 'LeftUpLeg', 'RightUpLeg',...
        'LeftLeg', 'RightLeg', 'LeftFoot', 'RightFoot'};
    m = containers.Map(key, val);
    for i=1:length(key)
        dataRaw.(key{i}) = dataRaw.(key{i})(1:nSamples,:);
        dataBVH.(val{i}) = dataBVH.(val{i})(1:nSamples,:)/1000;
    end
    dataBVH.posUnit = 'm';
    sIdx = 3; eIdx = length(dataBVH.Hips(:,1)) - 1;
    idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
    
    %% Calculate Orientation  
    qOri = {};
    qOri.x.PELV = dataBVH.qHips(sIdx:eIdx, :);
    qOri.x.LTIB = dataBVH.qLeftLeg(sIdx:eIdx, :);
    qOri.x.RTIB = dataBVH.qRightLeg(sIdx:eIdx, :);
    
    %% Position, Velocity, Acceleration
    gfrAcc = {};
    
    % gfrAcc from xsens
    gfrAcc.raw = {};
    gfrAcc.raw.MP = quatrotate(quatconj(dataRaw.Pelvis.ori), ...
                            dataRaw.Pelvis.acc) - [0 0 9.81];
    gfrAcc.raw.MP = gfrAcc.raw.MP(sIdx:eIdx,:);
    gfrAcc.raw.LA = quatrotate(quatconj(dataRaw.L_LowLeg.ori), ...
                            dataRaw.L_LowLeg.acc) - [0 0 9.81];
    gfrAcc.raw.LA = gfrAcc.raw.LA(sIdx:eIdx,:);
    gfrAcc.raw.RA = quatrotate(quatconj(dataRaw.R_LowLeg.ori), ...
                            dataRaw.R_LowLeg.acc) - [0 0 9.81];
    gfrAcc.raw.RA = gfrAcc.raw.RA(sIdx:eIdx,:);
    
    MIDPEL_act = [mean([dataBVH.LeftUpLeg(:,1) dataBVH.RightUpLeg(:,1)], 2),...
                  mean([dataBVH.LeftUpLeg(:,2) dataBVH.RightUpLeg(:,2)], 2),...
                  mean([dataBVH.LeftUpLeg(:,3) dataBVH.RightUpLeg(:,3)], 2)];
    gfr_vel_MP_act = diff(MIDPEL_act, 1, 1)*fs;
    gfr_vel_MP_act = [gfr_vel_MP_act(1,:); gfr_vel_MP_act];
    gfr_vel_LA_act = diff(dataBVH.LeftFoot, 1, 1)*fs;
    gfr_vel_LA_act = [gfr_vel_LA_act(1,:); gfr_vel_LA_act];
    gfr_vel_RA_act = diff(dataBVH.RightFoot, 1, 1)*fs;
    gfr_vel_RA_act = [gfr_vel_RA_act(1,:); gfr_vel_RA_act];
           
    x0_pos_MP = MIDPEL_act(sIdx,:);
    x0_pos_LA = dataBVH.LeftFoot(sIdx,:);
    x0_pos_RA = dataBVH.RightFoot(sIdx,:);
    x0_vel_MP = gfr_vel_MP_act(sIdx,:);
    x0_vel_LA = gfr_vel_LA_act(sIdx,:);
    x0_vel_RA = gfr_vel_RA_act(sIdx,:);
    
    gfrAcc.sim = {};
    gfrAcc.sim.MP = [0 0 0; diff(MIDPEL_act, 2, 1)*fs*fs];
    gfrAcc.sim.MP = gfrAcc.sim.MP(sIdx:eIdx,:);
    gfrAcc.sim.LA = [0 0 0; diff(dataBVH.LeftFoot, 2, 1)*fs*fs];
    gfrAcc.sim.LA = gfrAcc.sim.LA(sIdx:eIdx,:);
    gfrAcc.sim.RA = [0 0 0; diff(dataBVH.RightFoot, 2, 1)*fs*fs];
    gfrAcc.sim.RA = gfrAcc.sim.RA(sIdx:eIdx,:);
    
    %% UWB measurements
    %  Simulate uwb measurement by generating pairwise combinations, using the
    %  origin of each bone segment as the root point
    uwb_mea = struct;
    
    uwb_mea.left_tibia_mid_pelvis = vecnorm((MIDPEL_act-dataBVH.LeftFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.mid_pelvis_right_tibia = vecnorm((MIDPEL_act-dataBVH.RightFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);
    uwb_mea.left_tibia_right_tibia = vecnorm((dataBVH.RightFoot-dataBVH.LeftFoot), 2, 2) ...
        + normrnd(0, 0.02, [nSamples, 1]);

    d_pelvis = norm(dataBVH.RightUpLeg(sIdx,:) - dataBVH.LeftUpLeg(sIdx,:));
    d_rfemur = norm(dataBVH.RightUpLeg(sIdx,:) - dataBVH.RightLeg(sIdx,:));
    d_lfemur = norm(dataBVH.LeftUpLeg(sIdx,:) - dataBVH.LeftLeg(sIdx,:));
    d_rtibia = norm(dataBVH.RightLeg(sIdx,:) - dataBVH.RightFoot(sIdx,:));
    d_ltibia = norm(dataBVH.LeftLeg(sIdx,:) - dataBVH.LeftFoot(sIdx,:));
    
    %% Run Experiment
    actBody = dataBVH.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', ...
                         'xyzColor', {'m', 'y', 'c'}}); 
    actBodyRel = actBody.changeRefFrame('MIDPEL');
            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        csGfrAcc = gfrAcc.(cs.accData);
        csQOri = qOri.(cs.oriData);
        
        % step detection
        VAR_WIN  = floor(fs*cs.stepDetectWindow); % NUM_SAMPLES
        ACC_VAR_THRESH = cs.stepDetectThreshold;

        movVarAcc_pelvis = movingvar(sqrt( sum(csGfrAcc.MP .^2, 2)), VAR_WIN);
        bIsStatMP = movVarAcc_pelvis < 0;
        movVarAcc_lankle = movingvar(sqrt( sum(csGfrAcc.LA .^2, 2)), VAR_WIN);
        bIsStatLA = movVarAcc_lankle < ACC_VAR_THRESH;
        movVarAcc_rankle = movingvar(sqrt( sum(csGfrAcc.RA .^2, 2)), VAR_WIN);
        bIsStatRA = movVarAcc_rankle < ACC_VAR_THRESH;
        
%         try
            if cs.est == 'ekfv3'
                x0 = [x0_pos_MP x0_vel_MP zeros(1,4) ...
                      x0_pos_LA x0_vel_LA zeros(1,4) ...
                      x0_pos_RA x0_vel_RA zeros(1,4)]';
                v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigmaQAccMP', cs.sigmaQAcc, ...
                    'sigmaQAccLA', cs.sigmaQAcc, 'sigmaQAccRA', cs.sigmaQAcc);
                
                [ x_pri_v2, x_pos_v2, t_dat_v2 ] = grlib.est.kf_3_kmus_v3( ...
                    x0, cs.P, csGfrAcc.MP, bIsStatMP, csQOri.PELV, ...
                    csGfrAcc.LA, bIsStatLA, csQOri.LTIB, ...
                    csGfrAcc.RA, bIsStatRA, csQOri.RTIB, ...
                    d_pelvis, d_lfemur, d_rfemur, d_ltibia, d_rtibia, ...
                    uwb_mea, v3Options);
                
                estBody = grlib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
                   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'vicon', ...
                   'xyzColor', {'r', 'g', 'b'}, ...
                   'MIDPEL', x_pos_v2(idx0, 1:3), ...
                   'LFEP', t_dat_v2.LFEP(idx0, :), ...
                   'LFEO', t_dat_v2.LFEO(idx0, :), ...
                   'LTIO', x_pos_v2(idx0, 11:13), ...
                   'RFEP', t_dat_v2.RFEP(idx0, :), ...
                   'RFEO', t_dat_v2.RFEO(idx0, :), ...
                   'RTIO', x_pos_v2(idx0, 21:23), ...
                   'qRPV', x_pos_v2(idx0, 7:10), ...
                   'qLTH', t_dat_v2.qLTH(idx0, :), ...
                   'qRTH', t_dat_v2.qRTH(idx0, :), ...
                   'qLSK', x_pos_v2(idx0, 17:20), ...
                   'qRSK', x_pos_v2(idx0, 27:30));
               
                estState = x_pos_v2;
                estState2 = t_dat_v2;
                actState = [MIDPEL_act(sIdx:eIdx,:) ...
                   gfr_vel_MP_act(sIdx:eIdx,:) qOri.x.PELV...
                   dataBVH.LeftFoot(sIdx:eIdx,:) ...
                   gfr_vel_LA_act(sIdx:eIdx,:) qOri.x.LTIB...
                   dataBVH.RightFoot(sIdx:eIdx,:) ...
                   gfr_vel_RA_act(sIdx:eIdx,:) qOri.x.RTIB];
                estAcc = [csGfrAcc.MP csGfrAcc.LA csGfrAcc.RA];
                actAcc = [gfrAcc.sim.MP gfrAcc.sim.LA gfrAcc.sim.RA];
            end
            
            estBodyRel = estBody.changeRefFrame('MIDPEL');
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'actBody', 'estState', 'actState', ...
                     'estState2', 'estAcc', 'actAcc')
            end
    %         results(resultsIdx) = estBody.diffRMSE(actBody);
            results0 = estBodyRel.diffRMSE(actBodyRel);
%         catch
%             results0 = actBody.diffRMSE(nan);
%         end
        
        results0.name = name;
        results0.label = cs.label;
        results0.runtime = cputime-t0;
        results(resultsIdx) = results0;
        display(sprintf("Index %3d/%3d: Running time: %.4f", resultsIdx, setupN, cputime-t0));
        resultsIdx = resultsIdx + 1;
    end
end