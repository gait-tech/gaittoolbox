function varargout = quaternion2rotationmatrix(varargin)
%QUATERNION2ROTATIONMATRIX Converts a quaternion orientation to a rotation matrix
%
%   R = quaternion2rotationmatrix(q)
%
%   see: 'Representing Attitude Euler Angles, Unit Quaternions, and Rotation Vectors'
%
%   | q0^2+q1^2-q2^2-q3^2      2q1q2+2q0q3              2q1q3-2q0q2      |
%   |  
%   |     2q1q2-2q0q3      q0^2-q1^2+q2^2-q3^2          2q2q3+2q0q1      |
%   |
%   |     2q1q3+2q0q2          2q2q3-2q0q1          q0^2-q1^2-q2^2+q3^2  |      
%
% disp('quaternion2rotationmatrix:= see Representing Attitude Euler Angles')
%
% Rotation matrix may also be thought of as the matrix of basis vectors
% that define the world and body-fixed coordinates expressed in world
% coordinates. 
% COLUMNS basis vectors of world coordinates in body-fixed coordinates
% ROWS    basis vectors of body-fixed coordinates expressed in world frame
bRows = true;
if nargin == 1
    q = varargin{1};
elseif nargin == 2;
    q = varargin{1};
    outputType = varargin{2};
    % 'row' or 'column'
    if strcmp(outputType,'column') 
        bRows = false;
    end
else 
    error('Not supported');
end
    % row 1
    R(1,1,:) = q(:,1).^2 + q(:,2).^2 - q(:,3).^2 - q(:,4).^2;
    R(1,2,:) = 2.*(q(:,2).*q(:,3) + q(:,1).*q(:,4));
    R(1,3,:) = 2.*(q(:,2).*q(:,4) - q(:,1).*q(:,3) );
    % row 2
    R(2,1,:) = 2.*(q(:,2).*q(:,3) - q(:,1).*q(:,4));
    R(2,2,:) = q(:,1).^2 - q(:,2).^2 + q(:,3).^2 - q(:,4).^2;
    R(2,3,:) = 2.*(q(:,3).*q(:,4) + q(:,1).*q(:,2));
    % row 3
    R(3,1,:) = 2.*(q(:,2).*q(:,4) + q(:,1).*q(:,3) );
    R(3,2,:) = 2.*(q(:,3).*q(:,4) - q(:,1).*q(:,2));
    R(3,3,:) = q(:,1).^2 - q(:,2).^2 - q(:,3).^2 + q(:,4).^2;
    
if nargout == 1
    varargout{1} = R;
elseif nargout == 3
    if bRows
%         fprintf('Basis vectors of body coordinates in world frame\n');
        xhat = R(1,:,:);
        yhat = R(2,:,:);
        zhat = R(3,:,:);
        varargout{1} = reshape(xhat,3,[])';
        varargout{2} = reshape(yhat,3,[])';
        varargout{3} = reshape(zhat,3,[])';
    else
%         fprintf('Basis vectors of world frame in body-fixed coordinates\n');        
        xhat = R(:,1,:);
        yhat = R(:,2,:);
        zhat = R(:,3,:);        
        varargout{1} = reshape(xhat,3,[])';
        varargout{2} = reshape(yhat,3,[])';
        varargout{3} = reshape(zhat,3,[])';
    end
else
    error('Not Supported');
end

end

