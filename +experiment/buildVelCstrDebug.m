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
function sensors = buildVelCstrDebug(sensors, body, suffix)
    if nargin <= 2
        suffix = '';
    end
    
    vel = body.calcJointVel({'LFEO', 'RFEO', 'LFEP', 'RFEP'});
    avel = body.calcSegAngVel({'qLSK', 'qRSK'}, 'W');
    
    sensors.(sprintf('LT_v_dLK%s', suffix)) = quatrotate(body.qLTH, vel.LFEP-vel.LFEO-cross(avel.qLSK, body.LFEP-body.LFEO));
    sensors.(sprintf('RT_v_dRK%s', suffix)) = quatrotate(body.qRTH, vel.RFEP-vel.RFEO-cross(avel.qRSK, body.RFEP-body.RFEO));
end