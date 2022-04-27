function [control] = NMPC_soft(h, v, V, hd, k, vc, args, solver ,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
H = [h;v;V];
args.p(1:9) = H;

for i = 1:N
    args.p(9*i+1:9*i+9)=hd(:,k+i);
end

args.x0 = reshape(vc',4*N,1);
sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
    'lbg', args.lbg, 'ubg', args.ubg, 'p', args.p);
control = reshape(full(sol.x)',4,N)';

end

