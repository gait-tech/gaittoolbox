function [ b ] = bernoullinumber( k )
% BERNOULLINUMBER Generate the kth bernoulli number
% From: http://www.mathworks.com/matlabcentral/fileexchange/7257-bernoulli-numbers

if k == 0 
    b = 1;
elseif k == 1
    b = -1/2;
elseif mod(k,2)~= 0 
    b = 0;
else
    c = 1./factorial(2:k+1);
    b = (-1)^k .*factorial(k).* ...
    det(toeplitz(c, [c(1) 1 zeros(1, k-2)]));
end

end

