function [ul_ref, um_ref, un_ref, w_ref] = ref_velocities(hxdp, hydp, hzdp, psidp, psid)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for k =1:length(hxdp)
    J = [cos(psid(k)) -sin(psid(k))    0    0;...
        sin(psid(k))  cos(psid(k))    0    0;...
        0        0        1    0;...
        0        0        0    1];
    hp = [hxdp(k);hydp(k);hzdp(k);psidp(k)];
    
    uref = inv(J)*hp;
    
    ul_ref(k)= uref(1);
    um_ref(k)= uref(2);
    un_ref(k)= uref(3);
    w_ref(k)= uref(4);
end

end

