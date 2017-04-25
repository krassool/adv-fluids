function [u,v] = flow_field_cyl_1_0(Xj,Yj,q,xp,yp)

% make a mesh of x and y positions

Xmj = 0.5*(Xj(1) + Xj(2)); % midpoints
Ymj = 0.5*(Yj(1) + Yj(2));
% equation (13)
rij = sqrt((Xmj - xp).^2 + (Ymj - yp).^2);
% equation (14) and (15)
Phi = atan2((Yj(2) - Yj(1)),(Xj(2) - Xj(1)));
beta = atan2((yp - Ymj),(xp - Xmj));
omega = beta - Phi;
% equations (16) and (17)
x0p = rij.*cos(omega);
y0p = rij.*sin(omega);
% equation (11 & 12)
S = sqrt((Xj(2) - Xj(1)).^2 + (Yj(2) - Yj(1)).^2);
a = -S/2;
b = S/2;
% equations (18) and (19)
vprime = (q./(2*pi)).*(atan(((S/2) - x0p)./y0p) - atan((-(S/2) - x0p)./y0p));
uprime = (q./(2*pi)).*((-log((y0p.^2 + ((S.^2)/4) - (S.*x0p) + x0p.^2))./2) ...
    - (-log((y0p.^2 + ((S.^2)/4) + (S.*x0p) + x0p.^2))./2));
% equations (21) and (22)
v = vprime.*cos(Phi) + uprime.*sin(Phi);
u = uprime.*cos(Phi) - vprime.*sin(Phi);

end
