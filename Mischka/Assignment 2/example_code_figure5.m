% this is an example code to produce the right-hand-plot of figure 5

% make a mesh of x and y positions
[xp, yp] = meshgrid([-10.1:.5:10.1], [-10.1:.5:10.1]);
q = 1; % source panel strength density

Xj = [-1 -4]; % start and end points x and y coords
Yj = [4 1];
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

% plot the panel
plot(Xj, Yj, 'b-', 'LineWidth', 2.5)
hold on
% plot the vectors
quiver(xp, yp, u, v, 4, 'k')
daspect([1 1 1])
axis([-10 10 -10 10])
xlabel('x')
ylabel('y')

figure
[u2, v2] = source_panel_field(q, Xj, Yj, xp, yp);
% plot the panel
plot(Xj, Yj, 'b-', 'LineWidth', 2.5)
hold on
% plot the vectors
quiver(xp, yp, u2, v2, 4, 'k')
daspect([1 1 1])
axis([-10 10 -10 10])
xlabel('x')
ylabel('y')

