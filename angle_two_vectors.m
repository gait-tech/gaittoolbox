function [ cos_theta,angle_rad ] = angle_two_vectors( u,v )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

[N_U,~] = size(u);
[N_V,~] = size(v);

if N_U ~= N_V
end
cos_theta = NaN(N_U,1);
angle_rad = NaN(N_U,1);
for i=1:N_U
    cos_theta(i)= dot(u(i,:),v(i,:))/(norm(u(i,:))*norm(v(i,:)));
    angle_rad(i)    = acos(cos_theta(i));
end



end

