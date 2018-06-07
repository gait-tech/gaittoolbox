% ======================================================================
%> @brief Returns the difference between angle A and angle B
%>
%> Returns min(|angle A - angle B|, |angle A + angle B|)
%>
%> @param A angle input 1 
%> @param B angle input 2 
%>
%> @retval ret difference between angle A and B
% ======================================================================
function ret = anglediff(A, B)
    validateattributes(A, {'numeric'}, {'2d'})
    validateattributes(B, {'numeric'}, {'2d'})
    
    sizeA = size(A);
    ret = zeros(sizeA);
    
    for i=1:sizeA(1)
        for j=1:sizeA(2)
            if abs(A(i,j)-B(i,j)) < abs(A(i,j)+B(i,j)) 
                ret(i,j) = A(i,j)-B(i,j);
            else
                ret(i,j) = A(i,j)+B(i,j);
            end
        end
    end
end