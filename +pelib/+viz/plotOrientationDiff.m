function plotOrientationDiff(body1, body2, parts)
	% Plot body orientation (body2 - body1)
	%
	% :param body1: Base Body instance
	% :param body2: Body instance(s) to be compared to
	% :param parts: String(s) of body point(s) to be plotted.
	%
	% .. Author: - Luke Sy (UNSW GSBME)

    if ~iscell(body2)
        body2 = {body2};
    end
    if ~iscell(parts)
        parts = {parts};
    end
    
    n = length(parts);
    
    for i=1:n
        actl = quat2eul(body1.(parts{i}));
        if body1.oriUnit == 'deg'
            actl = actl * 180 / pi;
        end
        t = 1:length(actl(:,1));
        
        subplot(n,1,i); hold on;
        
        for j=1:length(body2)
            pred = quat2eul(body2{j}.(parts{i}));
            if body1.oriUnit == 'deg'
                pred = pred * 180 / pi;
            end
            
            xErr = (actl(:,1)-pred(:,1));
            yErr = (actl(:,2)-pred(:,2));
            zErr = (actl(:,3)-pred(:,3));
            tErr = vecnorm([xErr yErr zErr], 2, 2);

            plot(t,xErr, strcat(body2{j}.xyzColor{1}, body2{j}.lnSymbol), ...
                 t,yErr, strcat(body2{j}.xyzColor{2}, body2{j}.lnSymbol), ...
                 t,zErr, strcat(body2{j}.xyzColor{3}, body2{j}.lnSymbol), ...
                 t,tErr, strcat('k', body2{j}.lnSymbol));
        end
        title(parts{i});
        xlabel('Time');
        ylabel(strcat('Error (', body2{1}.oriUnit, ')'));
        legend('x','y','z','total');
    end
end