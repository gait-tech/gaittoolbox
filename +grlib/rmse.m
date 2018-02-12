% ======================================================================
%> @brief Returns the rmse of measurements
%>
%> @param objs objects
%>
%> @retval out rmse of objs
% ======================================================================
function out = rmse(objs)
    out = sqrt(nanmean(objs.^2, 1));
end