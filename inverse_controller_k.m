function [control_law] = inverse_controller_k(h, hd_k, hd, L, ts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
he = hd - h;

%% GAIN MATRICES OF THE SYSTEM
K1 = 0.97*eye(size(h,1));

J = drone_jacobian(h, L);

control_law = pinv(J)*((hd_k-K1*he-h)/ts);
control_law =  control_law';
end

