% ======================================================================
%> @brief build state debug
%> @author Luke Sy (UNSW GSBME)
%> @date 7 Aug 2019
%>
%> @param state estimate state
%> @param state2 debug data of estimator
%> @param algo algorithm string code
%>
%> @retval newState expanded state
% ======================================================================
function newState = buildStateDebug(state, state2, algo)
    newState = false;
    if strcmp(algo, 'lieekfv1')
        newState = struct('W_T_PV', repelem(state.W_T_PV, 1, 1, 3), ...
                          'W_T_LS', repelem(state.W_T_LS, 1, 1, 3), ...
                          'W_T_RS', repelem(state.W_T_RS, 1, 1, 3), ...
                          'vec', repelem(state.vec, 1, 3));
        n = size(newState.vec, 2);
        
        for i = {'W_T_PV', 'W_T_LS', 'W_T_RS'}
            sname = i{1};
            newState.(sname)(:,:,1:3:n) = state2.xhatPri.(sname);
            newState.(sname)(:,:,2:3:n) = state2.xhatPos.(sname);
        end
        newState.vec(:,1:3:n) = state2.xhatPri.vec;
        newState.vec(:,2:3:n) = state2.xhatPos.vec;
    end
end