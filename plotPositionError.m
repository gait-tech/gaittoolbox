% ======================================================================
%> @brief Plot the error between the predicted and actual position.
%> @author Luke Sy
%>
%> Plot the error between the predicted and actual position.
%> 
%> @param pred predicted position
%> @param actl actual position
%> @param dim if 1 dimension equals (n x 3), else dimension equals (3 x n)
%>
%> @retval p plot object
% ======================================================================
function p = plotPositionError(pred, actl, dim)
    if nargin <= 2
        dim = 1;
    end
        
    validateattributes(pred, {'numeric'}, {'2d'});
    validateattributes(actl, {'numeric'}, {'2d'});
    
    if dim ~= 1
        pred = pred';
        actl = actl';
    end
    
    t = 1:length(pred(:,1));
    xErr = (actl(:,1)-pred(:,1));
    yErr = (actl(:,2)-pred(:,2));
    zErr = (actl(:,3)-pred(:,3));
    tErr = vecnorm([xErr yErr zErr], 2, 2);
    
    p = plot(t,xErr,'-r',...
             t,yErr,'-g',...
             t,zErr,'-b',...
             t,tErr,'-k');
    xlabel('Time'); ylabel('Error');
    legend('x','y','z','total');
end