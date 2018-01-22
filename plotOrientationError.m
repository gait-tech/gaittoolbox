% ======================================================================
%> @brief Plot the error between the predicted and actual orientation.
%> @author Luke Sy
%>
%> Plot the error between the predicted and actual orientation.
%> Took the difference between the euler angle representation ('ZYX') of 
%> the 2 orientations 
%> 
%> @param pred predicted position
%> @param actl actual position
%> @param dim if 1 dimension equals (n x 3), else dimension equals (3 x n)
%>
%> @retval p plot object
% ======================================================================
function p = plotOrientationError(pred, actl, dim)
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
    predEul = quat2eul(pred);
    actlEul = quat2eul(actl);
    
    xErr = anglediff(predEul(:,1), actlEul(:,1))*180/pi;
    yErr = anglediff(predEul(:,2), actlEul(:,2))*180/pi;
    zErr = anglediff(predEul(:,3), actlEul(:,3))*180/pi;
    tErr = vecnorm([xErr yErr zErr], 2, 2);
    
    p = plot(t,xErr,'-r',...
             t,yErr,'-g',...
             t,zErr,'-b',...
             t,tErr,'-k');
    xlabel('Time'); ylabel('Error (Degree)');
    legend('x','y','z','total');
end