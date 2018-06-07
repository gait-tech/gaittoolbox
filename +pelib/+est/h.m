function [hp] = h(x,fs)
nDim = 3;
nSense = 3;
nStpSns = 9;
nSt = length(x);
hp = zeros(nDim*nSense,1);
k = 1;
for i = 0:nStpSns:nSt-nDim
    for j = 1:nSense
        hp(k) = x(i+j+2*nDim);
        k = k+1;
    end
end
%add noise
%hp = hp + diag(R);