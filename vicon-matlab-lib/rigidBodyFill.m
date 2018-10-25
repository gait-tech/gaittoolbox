% ======================================================================
%> @brief Rigid body fill
%> @author Luke Sy (UNSW GSBME)
%> @date 9 Oct 2018
%>
%> Training
%> base1, base2, base3 = input points to generate rigid body model
%> baseTarget = output point of rigid body model
%> 
%> Actual
%> est1, est2, est3 = input point to generate rigid body model
%> estTarget = output point of rigid body model
%>
%> @param base1 1 x 3 base vector (origin)
%> @param base2 1 x 3 base vector
%> @param base3 1 x 3 base vector
%> @param baseTarget 1 x 3 target base vector
%> @param est1 1 x 3 vector (origin)
%> @param est2 1 x 3 vector
%> @param est3 1 x 3 vector
%>
%> @retval estTarget 1 x 3 target vector
% ======================================================================
function estTarget = rigidBodyFill(base1, base2, base3, baseTarget, ...
                                   est1, est2, est3)
    baseR = getRotationMatrix(base1, base2, base3);
    v_GCS = (baseTarget - base1)'; % 3 x 1 target vector
    v_LCS = baseR'*v_GCS; % 3 x 1 vector. components of basis X, Y, Z
                          % v in local coordinate system
               
    n = size(est1, 1);
    estTarget = zeros(n, 3);
    for i=1:n
        estR = getRotationMatrix(est1(i,:), est2(i,:), est3(i,:));
        estTarget(i,:) = (estR*v_LCS)' + est1(i,:);
    end
end