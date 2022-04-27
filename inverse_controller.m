function [control_law] = inverse_controller(h, hd, hdp, k1, k2, L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
he = hd - h;

%% GAIN MATRICES OF THE SYSTEM
K1 = k1*eye(size(h,1));
K2 = k2*eye(size(h,1));

J = drone_jacobian(h, L);

control_law = pinv(J)*(hdp+K2*tanh(pinv(K2)*K1*he));
control_law =  control_law';
end

