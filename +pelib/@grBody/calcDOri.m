function out = calcDOri(obj, ref, targetSeg)
	% Calculate d_ori of lower body grBody as defined in the Marcard paper
	% 
	% T. von Marcard, B. Rosenhahn, M. J. Black, G. P.-M. (2017). 
	% Sparse Inertial Poser: Automatic 3D Human Pose Estimation from Sparse IMUs. 
	% Eurographics, 36(2)
	%
	% :param obj: this grBody
	% :param ref: reference grBody to be compared with
    % :param targetSeg: segments to be computed (usually occluded)
	%
	% :return: out - array of d_ori with respect to time
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 12/07/18

    if nargin <= 2, targetSeg = {'qLTH', 'qRTH'}; end
    rN = size(obj.qLTH, 1); cN = length(targetSeg);
    dori = zeros(rN, cN);
    
    cIdx = 0;
    for i=1:cN
        n = targetSeg{i};
        if (isempty(obj.(n)) || isempty(ref.(n))), continue; end
            
        buf = calcEul(quat2rotm(quatmultiply(obj.(n), quatconj(ref.(n)))));
        cIdx = cIdx + 1;
        dori(:, cIdx) = vecnorm(rad2deg(buf), 2, 2);
    end
    out = mean(dori(:,1:cIdx), 2);
end

function eul = calcEul(R)
    n = size(R, 3);
    eul = zeros(n, 3);
    for i=1:n
        if any(any(isnan(R(:,:,i))))
            eul(i, :) = [nan nan nan];
        else
            R2 = logm(R(:,:,i));
            eul(i, :) = [R2(3,2) R2(2,1) R2(1,3)];
        end
    end
end