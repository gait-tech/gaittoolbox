function [ rotatedVecs ] = rotateVecsByQuats( vecs,quats )
%ROTATEVECSBYQUATS rotate a matrix of vectors by a corresponding matrix of
%quaternions
%   [ rotatedVecs ] = rotateVecsByQuats( vecs,quats )
[numVecs, vecDim]  = size(vecs);
[numQuats,quatDim] = size(quats);
rotatedVecs = zeros(numQuats,3);
if quatDim ~= 4, error('Quaternions must have 4 dimensions');end
if vecDim ~=3, error('Vectors must be 3 dimensional');end

if numVecs ~= numQuats,
    if numVecs == 1
        % rotate the vector by each quaternion    
        fprintf('Rotating the same Vector by each Quaternion\n');
        for i = 1:numQuats
            rotatedVecs(i,:) = quatrot(vecs,quats(i,:));
        end
    end
    if numQuats == 1
        rotatedVecs = zeros(numVecs,3);
        fprintf('Rotating each Vector by the same Quaternion\n');
        for i = 1:numVecs
            rotatedVecs(i,:) = quatrot(vecs(i,:),quats);
        end        
    else
        fprintf('Unequal number of Quaternions and Vectors\n- rotating each vector by each quaternion\n');
        rotatedVecs = zeros(numQuats*numVecs,3);
        for i = 1:numQuats
            for j=1:numVecs
                rotatedVecs( (i-1)*numVecs +j ,:) = quatrot(vecs(j,:),quats(i,:));
            end
        end
    end
else
    % rotate each vector by its matching quaternion    
    for i=1:numVecs
        rotatedVecs(i,:) = quatrot(vecs(i,:),quats(i,:));
    end    
end
end

function rVec = quatrot(v,q)

% Use quaternion to rotate vector
% 
% v: vector to rotate
% q: unit quaternion to perform rotation around vector (x,y,z) using 
% the quaternion q = w + xi +yj +zk
%
% Ensure column vectors
% q = q(:);
% v = v(:);

% Check size
if length(q)~=4
    error('Quaternion must contain four elements.');
end

% Ensure q is unit quaternion
normq = 1/realsqrt(sum(q.^2));
q = q.*normq;

% vrot = qvq*
vrot = quatmultiply(quatmultiply(q,[0  v]),quatinv(q));
rVec = vrot(2:4);

end
