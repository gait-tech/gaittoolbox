% ======================================================================
%> @brief Calculate distal segment rotation matric from proximal segment and joint angle
%>
%> @param prox quaternion orientation (n x 4) of the proximal segment
%> @param theta joint angle in X, Y, Z axis
%> @param seq (default: YX'Z'')
%>
%> @retval theta [x y z] radians 
% ======================================================================
function theta = calcJointAngles(prox, dist, seq)
    if nargin <= 2
        seq = 'YXZ';
    end
    
%     pCS = quat2rotm(prox);
%     dCS = quat2rotm(dist);
    relori = quatmultiply(quatconj(prox), dist);
%     relori = quatmultiply(quatconj(dist), prox);
%     [theta_y theta_x theta_z] = quat2angle(relori, seq);
    % Cardan angles calculation
%     theta_x = reshape(-asin(dot(dCS(:,3,:), pCS(:,2,:))), [], 1);
%     theta_y = reshape(asin(dot(dCS(:,3,:), pCS(:,1,:))), [], 1)./cos(theta_x);
%     theta_z = reshape(asin(dot(dCS(:,1,:), pCS(:,2,:))), [], 1)./cos(theta_x);

    [theta1, theta2, theta3] = quat2angle(relori, seq);
    theta = [theta1 theta2 theta3];
%     theta = [theta_x theta_y theta_z];
end