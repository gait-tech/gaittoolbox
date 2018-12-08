% ======================================================================
%> @brief Calculate d_ori of lower body grBody as defined in the Marcard
%> paper
%> @author Luke Sy (UNSW GSBME)
%> @date 7 Dec 2018
%>
%> T. von Marcard, B. Rosenhahn, M. J. Black, G. P.-M. (2017). 
%> Sparse Inertial Poser: Automatic 3D Human Pose Estimation from Sparse IMUs. 
%> Eurographics, 36(2)
%>
%> @param obj this grBody
%> @param ref reference grBody to be compared with
%>
%> @retval out array of d_ori
% ======================================================================
function out = calcDOri(obj, ref)
    nameList = {'qLTH', 'qRTH'};
    doriList = {};
    
    for i=1:length(nameList)
        n = nameList{i};
        doriList.(n) = cellfun(@(x) logm(x), ...
            num2cell(quat2rotm(quatmultiply(obj.qLTH, quatconj(ref.qLTH))), [1 2]), ...
            'UniformOutput', false)
    end
    out = [];
end