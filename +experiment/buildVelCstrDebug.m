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
        LVcstrByVel =  quatrotate(body.qLTH, LHIPVel' - LKNEVel' - cross(state.vec(13:15, :)', (body.LFEP-body.LFEO)));
        RVcstrByVel =  quatrotate(body.qRTH, RHIPVel' - RKNEVel' - cross(state.vec(16:18, :)', (body.RFEP-body.RFEO)));
        
        sensors.(sprintf('LHIPVel%s', suffix)) = LHIPVel';
        sensors.(sprintf('RHIPVel%s', suffix)) = RHIPVel';
        sensors.(sprintf('LKNEVel%s', suffix)) = LKNEVel';
        sensors.(sprintf('RKNEVel%s', suffix)) = RKNEVel';
        sensors.(sprintf('LVcstrByVel%s', suffix)) = LVcstrByVel;
        sensors.(sprintf('RVcstrByVel%s', suffix)) = RVcstrByVel;
        
        % == start ==
        % copy paste parts of e23bVelCstrValidation.m here
        % == end ==
    end
end