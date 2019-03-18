% ======================================================================
%> @brief Calculate the (relative) joint angle between 2 body segments
%>
%>
%> @param prox quaternion orientation (n x 4) of the proximal segment
%> @param angles joint angles in seq order
%> @param seq (default: YX'Z'')
%>
%> @retval dist quaternion orientation (n x 4) of the distal segment
% ======================================================================
function dist = calcDistR(prox, angles, seq)
    if nargin <= 2
        seq = 'YXZ';
    end
    
    relori = angle2quat(angles(:,1), angles(:,2), angles(:,3), seq);
    dist = quatmultiply(prox, relori);
end