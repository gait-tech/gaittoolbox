% ======================================================================
%> @brief generate additional information for debugging the velocity constraint
%> @author Luke Sy (UNSW GSBME)
%> @date 1 Aug 2019
%>
%>
%> @param sensors struct containing sensor data
%> @param body grBody input
%>
%> @retval out grBody
% ======================================================================
function sensors = buildVelCstrDebug(sensors, body, state, algo, suffix)
    if nargin <= 4
        suffix = '';
    end
    
    vel = body.calcJointVel({'LFEO', 'RFEO', 'LFEP', 'RFEP'});
    avel = body.calcSegAngVel({'qLSK', 'qRSK'}, 'W');

    %% calculate difference from position states
    sensors.(sprintf('LVcstrByPos%s', suffix)) = quatrotate(body.qLTH, vel.LFEP-vel.LFEO-cross(avel.qLSK, body.LFEP-body.LFEO));
    sensors.(sprintf('RVcstrByPos%s', suffix)) = quatrotate(body.qRTH, vel.RFEP-vel.RFEO-cross(avel.qRSK, body.RFEP-body.RFEO));
        
    if strcmp(algo, 'lieekfv1')
        %% calculate difference from vel and angvel states
        LHIPVel = state.vec(1:3, :) + cross(state.vec(10:12, :), (body.LFEP-body.MIDPEL)');
        RHIPVel = state.vec(1:3, :) + cross(state.vec(10:12, :), (body.RFEP-body.MIDPEL)');
        LKNEVel = state.vec(4:6, :) + cross(state.vec(13:15, :), (body.LFEO-body.LTIO)');
        RKNEVel = state.vec(7:9, :) + cross(state.vec(16:18, :), (body.RFEO-body.RTIO)');
        LVcstrByVel =  LHIPVel - LKNEVel - cross(state.vec(13:15, :), (body.LFEP-body.LFEO)');
        RVcstrByVel =  RHIPVel - RKNEVel - cross(state.vec(16:18, :), (body.RFEP-body.RFEO)');
        
        sensors.(sprintf('LHIPVel%s', suffix)) = LHIPVel';
        sensors.(sprintf('RHIPVel%s', suffix)) = RHIPVel';
        sensors.(sprintf('LKNEVel%s', suffix)) = LKNEVel';
        sensors.(sprintf('RKNEVel%s', suffix)) = RKNEVel';
        sensors.(sprintf('LVcstrByVel%s', suffix)) = quatrotate(body.qLTH, LVcstrByVel');
        sensors.(sprintf('RVcstrByVel%s', suffix)) = quatrotate(body.qRTH, RVcstrByVel');

        %% equation testing
        n = body.nSamples;
        compL1 = zeros(n,18);
        compL2 = zeros(n,18);
        compR1 = zeros(n,18);
        compR2 = zeros(n,18);
        
        compL1(:,1:3) = quatrotate(body.qLTH, state.vec(1:3, :));
        compL1(:,4:6) = -quatrotate(body.qLTH, state.vec(4:6, :));
        compL1(:,7:9) = 0;
        compL1(:,10:12) = quatrotate(body.qLTH, cross(state.vec(10:12, :), (body.LFEP-body.MIDPEL)'));
        compL1(:,13:15) = quatrotate(body.qLTH, cross(state.vec(13:15, :), (body.LFEP-body.LTIO)'));
        compL1(:,16:18) = 0;
        
        addpath('liese3lib');
        
        VcstrByVel2 = zeros([2*size(LVcstrByVel,1) size(LVcstrByVel,2)]);
        
        for k=1:size(LVcstrByVel,2)
            body2 = struct();
            body2.PV_d = norm(body.LFEP(1,:)-body.RFEP(1,:));
            body2.LS_d = norm(body.LFEO(1,:)-body.LTIO(1,:));
            body2.RS_d = norm(body.RFEO(1,:)-body.RTIO(1,:));
            D0.H2CT = [eye(3,3); zeros(1,3)]';
            D0.PV_p_LH = [0 body2.PV_d/2 0 1]'; 
            D0.PV_p_RH = [0 -body2.PV_d/2 0 1]';
            D0.LS_p_LK = [0 0 body2.LS_d 1]';    
            D0.RS_p_RK = [0 0 body2.RS_d 1]';

            n_LT = D0.H2CT*(state.W_T_PV(:,:,k)*D0.PV_p_LH - ...
                            state.W_T_LS(:,:,k)*D0.LS_p_LK);
            n_RT = D0.H2CT*(state.W_T_PV(:,:,k)*D0.PV_p_RH - ...
                            state.W_T_RS(:,:,k)*D0.RS_p_RK);
            W_vy_PV = body2.PV_d/2*state.W_T_PV(1:3, 2, k);

            W_ry_LT = state.W_T_LS(1:3, 2, k);
            W_ry_RT = state.W_T_RS(1:3, 2, k);
            W_rx_LS = state.W_T_LS(1:3, 1, k);
            W_rx_RS = state.W_T_RS(1:3, 1, k);

%             velcstrY = [W_ry_LT', -W_ry_LT', zeros(1,3) ...
%                 (hat(W_vy_PV)*W_ry_LT)', ...
%                 (-hat(n_LT)*W_ry_LT+W_rx_LS)', zeros(1,3);
%                 W_ry_RT', zeros(1,3), -W_ry_RT' ...
%                 (hat(-W_vy_PV)*W_ry_RT)', ...
%                 zeros(1,3), (-hat(n_RT)*W_ry_RT+W_rx_RS)' ];
            velcstrY = [W_ry_LT', -W_ry_LT', zeros(1,3) ...
                (hat(W_vy_PV)*W_ry_LT)', ...
                (-hat(n_LT)*W_ry_LT+body2.LS_d*W_rx_LS)', zeros(1,3);
                W_ry_RT', zeros(1,3), -W_ry_RT' ...
                (hat(-W_vy_PV)*W_ry_RT)', ...
                zeros(1,3), (-hat(n_RT)*W_ry_RT+body2.RS_d*W_rx_RS)' ];
            
            for l=1:3:18
                compL2(k,l+1) = velcstrY(1,l:l+2)*state.vec(l:l+2,k);
                compR2(k,l+1) = velcstrY(2,l:l+2)*state.vec(l:l+2,k);
            end
            
            W_vz_LS = body2.LS_d*state.W_T_LS(1:3, 3, k);
            W_vz_RS = body2.RS_d*state.W_T_RS(1:3, 3, k);

            velcstrZ = [n_LT', -n_LT', zeros(1,3) ...
                (hat(W_vy_PV)*n_LT)', (-hat(W_vz_LS)*n_LT)', zeros(1,3);
                n_RT', zeros(1,3), -n_RT' ...
                (hat(-W_vy_PV)*n_RT)', zeros(1,3), (-hat(W_vz_RS)*n_RT)' ];

            VcstrByVel2(:,k) = [0; 0; velcstrY * state.vec(:,k); ...
                                velcstrZ * state.vec(:,k)];
        end
        clf; pelib.viz.plotXYZ(100, LVcstrByVel', VcstrByVel2([1 3 5],:)');
    end
end