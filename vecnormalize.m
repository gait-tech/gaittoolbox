function [ column_norm ] = vecnormalize( matrix_row_vectors )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
column_norm = sqrt(sum(matrix_row_vectors.^2,2));

end

