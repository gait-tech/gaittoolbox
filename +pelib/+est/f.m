function [xp] = f_addN(x,fs)
nDim = 3;
nSense = 3;
nStpSns = 9;
nSt = length(x);
xp = zeros(nSt,1);
for i = 0:nStpSns:nSt-nDim
    for j = 1:nSense
        xp(i+j) = x(i+j)+(x(i+j+nDim)*1/fs)+(0.5*(1/fs)^2*x(i+j+2*nDim));...
        xp(i+j+nDim) = x(i+j+nDim)+1/fs*x(i+j+2*nDim);...
        xp(i+j+2*nDim) = x(i+j+2*nDim);
%     disp('i = ')
%     disp(i)
%     disp('j = ')
%     disp(j)
    end
end
%add noise
%xp = xp + diag(Q);

%         xp(1) = x(1)+(x(4)*1/fs)+(0.5*(1/fs)^2*x(7));...
%         xp(2) = x(2)+(x(5)*1/fs)+(0.5*(1/fs)^2*x(8));...
%         xp(3) = x(3)+(x(6)*1/fs)+(0.5*(1/fs)^2*x(9));...
%         xp(4) = x(4)+1/fs*x(7);...         
%         xp(5) = x(5)+1/fs*x(8);...
%         xp(6) = x(6)+1/fs*x(9);...
%         xp(7) = x(7);
%         xp(8) = x(8);
%         xp(9) = x(9);


