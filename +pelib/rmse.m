% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
% @brief Returns the rmse of measurements
%
% :param objs objects
%
% :return: out rmse of objs
% .. Author: - Luke Sy (UNSW GSBME) - 9/22/18
function out = rmse(objs)
    out = sqrt(nanmean(objs.^2, 1));
end