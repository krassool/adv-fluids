function [u,v] = vor_panel_on_point_vel(Xj,Yj,gam,xp,yp)
% this is changed from the first half of the assignment. This function uses
% a vortex panel to determine the effect on the point velocities in space.

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
a = -S/2 ; b = S/2;

% equations (18) and (19)
uprime = gam/(2*pi).*( atan((x0p-b)./y0p) - atan((x0p-a)./y0p)) ;
vprime = -gam/(2*pi).*( log((x0p-b).^2+y0p.^2)/2 - log((x0p-a).^2+y0p.^2)/2 ) ;

% equations (21) and (22)
f=0;
v = vprime.*cos(Phi+f) + uprime.*sin(Phi+f);
u = uprime.*cos(Phi+f) - vprime.*sin(Phi+f);
end
