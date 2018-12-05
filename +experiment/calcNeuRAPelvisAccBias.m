% ======================================================================
%> @brief Calculate pelvis accelerometer bias from sitting trial
%> @author Luke Sy
%> 
%>
%> @param fnameV loaded mocapdb.ViconBody
%> @param calibV2W quaternion (1 x 4) transforming vicon frame to world frame
%> @param calibW2V mocapdb.XsensBody transforming each sensor's world frame
%>                 to vicon frame
%> @param fnameX loaded mocapdb.BVHBody 
%> @param fnameS loaded mocapdb.XsensBody
%> @param startFrame frame number at which the algorithm will start
%> @param endFrame frame number at which the algorithm will end
% ======================================================================
function bias = calcNeuRAPelvisAccBias(dataS, dataV, ...
                        calibV2W, calibW2V, dataX, ...
                        startFrame, endFrame)
    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody', 'numeric'}, {});
    validateattributes(calibV2W, {'numeric'}, {});
    validateattributes(calibW2V, {'mocapdb.XsensBody', 'numeric'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 5 || startFrame < 0, startFrame = 100; end
    if nargin <= 6 || endFrame < 0, endFrame = inf; end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    bias = {};
    
    %% Generate vicon based inputs in world frame
    if ~isempty(dataV) & ~isempty(calibV2W)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        W__dataV = dataV.getSubset(1:nSamples).toWorldFrame(calibV2W);
        W__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        
        sIdx = max(W__dataV.getStartIndex()+1, startFrame);
        eIdx = min(length(W__dataV.PELV(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        
        % gfrAcc from vicon
        W__viconBody = W__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});  
        acc = W__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
        accMPRef = acc.MIDPEL(sIdx:eIdx,:);
        
        % gfrAcc from sparse
        accMPEst = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        accMPEst = accMPEst(sIdx:eIdx,:);
%         qAve = quaternion(W__dataS.Pelvis.ori(sIdx:eIdx, :));
        bias.w__v = quatrotate(W__dataS.Pelvis.ori(sIdx, :), ...
                            mean(accMPEst, 1) - mean(accMPRef, 1));
    end
    
    if ~isempty(dataV) & ~isempty(calibW2V)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        V__dataV = dataV.getSubset(1:nSamples);
        V__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        
        sIdx = max(V__dataV.getStartIndex()+1, startFrame);
        eIdx = min(length(V__dataV.PELV(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        
        %% position, velocity, acceleration
        V__viconBody = V__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});  
        acc = V__viconBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
        accMPRef = acc.MIDPEL(sIdx:eIdx,:);

        % gfrAcc from sparse
        accMPEst = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        accMPEst = accMPEst(sIdx:eIdx,:);
        accMPEst = quatrotate(quatconj(calibW2V.Pelvis.ori), accMPEst);
        
        q = quatmultiply(W__dataS.Pelvis.ori(sIdx, :), calibW2V.Pelvis.ori);
        bias.v__v = quatrotate(q, mean(accMPEst, 1) - mean(accMPRef, 1));
    end
    
    if ~isempty(dataX)
        nSamples = min(dataX.nSamples, dataS.nSamples);
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        
        dataX = dataX.toWorldFrame(qXsensV2W);
        W__dataX = dataX.getSubset(1:nSamples);
        W__dataX.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples).toViconFrame(calibW2V);

        sIdx = startFrame;
        eIdx = min(length(W__dataX.Hips(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
       
        % gfrAcc from xsens
        W__xsensBody = W__dataX.togrBody(1:nSamples, {'name', 'xsens', 'oriUnit', 'deg', ...
                             'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                             'xyzColor', {'m', 'y', 'c'}}); 
        acc = W__xsensBody.calcJointAcc({'MIDPEL', 'LTIO', 'RTIO'});
        accMPRef = acc.MIDPEL(sIdx:eIdx,:);
        
        % gfrAcc from sparse
        accMPEst = quatrotate(quatconj(W__dataS.Pelvis.ori), ...
                                W__dataS.Pelvis.acc) - [0 0 9.81];
        accMPEst = accMPEst(sIdx:eIdx,:);
        
        bias.w__x = quatrotate(W__dataS.Pelvis.ori(sIdx, :), ...
                            mean(accMPEst, 1) - mean(accMPRef, 1));
    end
    
end