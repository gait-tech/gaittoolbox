function [ varargout ] = quaternion2nautical( varargin )
%QUATERNION2NAUTICAL Calculates euler angles from unit quaternion as per
%described in "Representing Attitude: Euler Angles, Unit Quaternions, and 
%Rotation Vectors" - Author: James Diebel, Stanford University 20 October 
%2006. This sequence is often referred to as cardan angles, tait-bryan
%angles or nautical angles. They are also referred to as euler angles in
%aeronautics which creates the confusion.
%Phi   <-> roll  <-> bank     (rotation about x)
%Theta <-> pitch <-> attitude (rotation about y)
%Psi   <-> yaw   <-> heading  (rotation about z)
%[phi_roll,theta_pitch,psi_yaw] = quaternion2nautical( q )
%{'\phi;','\theta';'\psi'}
%Body fixed x axis - Forward direction is along the body
%Body fixed y-axis is starboard    (to the right)
%Body fixed z-axis points downward (to the ground) 
%The home position [0,0,0] is flat level pointing along the world x-axis.
%i.e. north. 
%Whilst defined for A NED reference frame, also suitable for an 
%Up-North-West reference frame in which:
%Body fixed x axis - Forward direction is along the body
%Body fixed y-axis is port (to the left/west)
%Body fixed z-axis points upward (to the sky) 

%Algorithm as per described in section 5.6.3 Euler Angles <= Rotation 
%Matrix, page 11. 
%
%Remember! Rotation matrices can be thought of as the matrix of basis 
%vectors that define the world and body-fixed coordinate systems. The rows
%of the rotation matrix are the basis vectors of the body-fixed coordinates
%expressed in world coordinates, the columns are the basis vectors of the
%world coordinates expressed in the body fixed coordinates
%
%NB: q must be a 4 column vector
%   [ nauticalAngles ] = quaternion2nautical( q )
%   where nauticalAngles = [phi,theta,psi] <-> [roll,pitch,yaw]
bDegrees = false;
if nargin >=1
    q = varargin{1};
end
if nargin >=2
    bDegrees = true;
end
[~,c]=size(q);
if c~=4
    q=q'; % take transpose to see if incorrectly orientated
    [~,c]=size(q);
    if c~=4
        error('Input is not a quaternion matrix or vector');
    end
end

zero  = 1;  one   = 2;  two   = 3;  three = 4;

% definitions are row major
r23 = 2*( q(:,two).*q(:,three)+q(:,zero).*q(:,one));
r33 = q(:,three).^2-q(:,two).^2-q(:,one).^2+q(:,zero).^2;

r13 = 2*( q(:,one).*q(:,three) - q(:,zero).*q(:,two));

r12 = 2*( q(:,one).*q(:,two) + q(:,zero).*q(:,three));
r11 = q(:,one).^2+q(:,zero).^2-q(:,three).^2-q(:,two).^2;

phi_roll    = atan2(r23,r33); 
theta_pitch = -asin(r13);
psi_yaw     = atan2(r12,r11); 

if bDegrees
    phi_roll    = rad2deg(phi_roll);
    theta_pitch = rad2deg(theta_pitch);
    psi_yaw     = rad2deg(psi_yaw);
end

if nargout <= 1
    varargout{1} = [phi_roll,theta_pitch,psi_yaw];
elseif nargout == 3 
    varargout{1} = phi_roll;
    varargout{2} = theta_pitch;
    varargout{3} = psi_yaw;
else
    error('Number of output arguments not supported');
end

