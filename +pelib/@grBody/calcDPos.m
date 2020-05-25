<<<<<<< HEAD
function out = calcDPos(obj, ref)
=======
function out = calcDPos(obj, ref, includeRoot)
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
	% Calculate d_pos of lower body grBody as defined in the Marcard paper
	% 
	% T. von Marcard, B. Rosenhahn, M. J. Black, G. P.-M. (2017). 
	% Sparse Inertial Poser: Automatic 3D Human Pose Estimation from Sparse IMUs. 
	% Eurographics, 36(2)
	%
	% :param obj: this grBody
	% :param ref: reference grBody to be compared with
<<<<<<< HEAD
	%
=======
	% :param includeRoot: [boolean] if root/pelvis is included
    %
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
	% :return: out - array of d_pos with respect to time
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 12/08/18

<<<<<<< HEAD
    nameList = {'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', 'RTIO'};

    rN = size(obj.LFEP, 1); cN = length(nameList);
    dpos = zeros(rN, cN);
    
    for i=1:cN
        n = nameList{i};
        dpos(:, i) = vecnorm(obj.(n)-ref.(n), 2, 2);
    end
    out = mean(dpos, 2);
=======
    if nargin <= 2, includeRoot = true; end
%     nameList = {'LFEP', 'LFEO', 'LTIO', 'RFEP', 'RFEO', 'RTIO'};
    if includeRoot
        nameList = obj.posList;
    else
        idx = strcmp(obj.posList, 'MIDPEL');
        nameList = obj.posList(~idx);
    end
    rN = size(obj.LFEP, 1); cN = length(nameList);
    dpos = zeros(rN, cN);
    
    cIdx = 0;
    for i=1:cN
        n = nameList{i};
        if isempty(obj.(n)) || isempty(ref.(n)), continue; end

        cIdx = cIdx + 1;
        dpos(:, cIdx) = vecnorm(obj.(n)-ref.(n), 2, 2);
    end
    out = mean(dpos(:,1:cIdx), 2);
>>>>>>> 8860699ab93014d7c72b14f3600fe1b99132d583
end