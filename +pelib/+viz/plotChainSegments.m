function p = plotChainSegments(Ts, p)
    % Suppose Ts = [T_1, T_2, T_3, ....] 
    % where T_i is a 4 x 4 homogenous matrix
    %
    % This function draws line between the translation components of
    % T_1, T_1*T_2, T_1*T_2*T_3, ...
    % 
	% :param Ts: 4 x 4 x n array of homogenenous matrices
    % :param p: [Optional] plot object, if available, will only update the line instead of creating a new one
    %
	% :return: p - plot object
	%
	% .. Author: - Luke Sy (UNSW GSBME, 20200303)
    if nargin <= 1
        p = [];
    end
    n = size(Ts,3);
    pt = zeros(3,n);
    Tbuf = eye(4);
    for i=1:n
        Tbuf = Tbuf*Ts(:,:,i);
        pt(:,i) = Tbuf(1:3,4);
    end
    
    if isempty(p)
        p = line(pt(1,:),pt(2,:),pt(3,:));
    else
        p.XData = pt(1,:);
        p.YData = pt(2,:);
        p.ZData = pt(3,:);
    end
end