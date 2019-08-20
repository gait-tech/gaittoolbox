% ======================================================================
%> @brief general samples for scatter plot
%> @author Luke Sy
%>
%> 
%> @param state system mean
%> @param P system convariance matrix
%> @param N number of samples to be generated
%> @param idx index
%> @param algo algo
%> @retval out struct of PV, LS, RS of populated samples
% ======================================================================

function out = e15sample(state, P, N, idx, dt, algo)
    if strcmp(algo, 'lieekfv1')
        P2 = (P(1:18,1:18,idx)+P(1:18,1:18,idx)')/2;
        r = mvnrnd(zeros(1,18), P2, N)';
        out = struct('PV', zeros(N,3), 'LS', zeros(N,3), 'RS', zeros(N,3));
        for i=1:N
            T = state.W_T_PV(:,:,idx)*vec2tran(r( 1: 6,i)*dt);
            out.PV(i,:) = T(1:3,4)';
            T = state.W_T_LS(:,:,idx)*vec2tran(r( 7:12,i)*dt);
            out.LS(i,:) = T(1:3,4)';
            T = state.W_T_RS(:,:,idx)*vec2tran(r(13:18,i)*dt);
            out.RS(i,:) = T(1:3,4)';
        end
    elseif strcmp(algo, 'lieekfv2')
        P2 = (P(1:18,1:18,idx)+P(1:18,1:18,idx)')/2;
        r = mvnrnd(zeros(1,18), P2, N)';
        out = struct('PV', zeros(N,3), 'LS', zeros(N,3), 'RS', zeros(N,3));
        for i=1:N
%             out.PV(i,:) = (state.W_R_PV(:,:,idx)*vec2rot(r(1:3,i)*dt)+r(10:12,i))';
%             out.LS(i,:) = (state.W_R_LS(:,:,idx)*vec2rot(r(4:6,i)*dt)+r(13:15,i))';
%             out.RS(i,:) = (state.W_R_RS(:,:,idx)*vec2rot(r(7:9,i)*dt)+r(16:18,i))';
            out.PV(i,:) = (state.vec(1:3,idx)+r(10:12,i)*dt)';
            out.LS(i,:) = (state.vec(4:6,idx)+r(13:15,i)*dt)';
            out.RS(i,:) = (state.vec(7:9,idx)+r(16:18,i)*dt)';
        end
    end
end