function [out1, out2] = calcDOrinobias(obj, ref)
	% Calculate d_ori of lower body grBody as defined in the Marcard paper
	% 
	% T. von Marcard, B. Rosenhahn, M. J. Black, G. P.-M. (2017). 
	% Sparse Inertial Poser: Automatic 3D Human Pose Estimation from Sparse IMUs. 
	% Eurographics, 36(2)
	%
	% :param obj: this grBody
	% :param ref: reference grBody to be compared with
	%
	% :return: out1 - array of d_ori no bias with respect to time
	% :return: out2 - array of bias
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 12/07/18

    nameList = {'qLTH', 'qRTH'};
    doriList = {};
    rN = size(obj.qLTH, 1); cN = length(nameList);
    dori = zeros(rN, cN);
    bias = zeros(cN, 3);
    for i=1:cN
        n = nameList{i};
        buf = calcEul(quat2rotm(quatmultiply(obj.(n), quatconj(ref.(n)))));
        bias(i, :) = nanmean(buf, 1);
        dori(:, i) = vecnorm(rad2deg(buf-bias(i, :)), 2, 2);
    end
    out1 = mean(dori, 2);
    out2 = mean(bias, 1);
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