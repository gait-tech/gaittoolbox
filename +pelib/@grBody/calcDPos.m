% ======================================================================
%> @brief Calculate d_pos of lower body grBody as defined in the Marcard
%> paper
%> @author Luke Sy (UNSW GSBME)
%> @date 8 Dec 2018
%>
%> T. von Marcard, B. Rosenhahn, M. J. Black, G. P.-M. (2017). 
%> Sparse Inertial Poser: Automatic 3D Human Pose Estimation from Sparse IMUs. 
%> Eurographics, 36(2)
%>
%> @param obj this grBody
%> @param ref reference grBody to be compared with
%>
%> @retval out array of d_pos with respect to time
% ======================================================================
function out = calcDPos(obj, ref)
    nameList = {'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', 'RTIO'};

    rN = size(obj.LFEP, 1); cN = length(nameList);
    dpos = zeros(rN, cN);
    
    for i=1:cN
        n = nameList{i};
        dpos(:, i) = vecnorm(obj.(n)-ref.(n), 2, 2);
    end
    out = mean(dpos, 2);
end