% MCEN90018: Advanced Fluid Dynamics - Assignment 2
% ------------------------------------------------------------------------
% Mischka Kamener  539030                           Last modified: 28/4/16
%
% Function based of "example_code_figure5.m" written by Nick Hutchins.
% Calculates the u and v velcoity at locations (xp, yp) due to a vortex
% panel with circulation gamma placed between (x(1) y(1)) and (x(2) y(2)).
function [u, v] = vortex_panel_field(gamma, x, y, xp, yp)

% Calculate midpoint
xm = [mean(x), mean(y)];

% Calculate distances
r = sqrt((xp - xm(1)).^2 + (yp - xm(2)).^2);
S = sqrt(diff(x).^2 + diff(y).^2);

% Calculate angles
phi   = atan2(diff(y), diff(x));
beta  = atan2(yp - xm(2), xp - xm(1));
omega = beta - phi;

% Calculate distances x0' and y0'
x0_d = r.*cos(omega);   y0_d = r.*sin(omega);

% Velocities in rotated coordinates u' and v'
v_d = (gamma/(2*pi))*(-0.5*log(y0_d.^2 + S^2/4 - x0_d.*S + x0_d.^2) + ...
    0.5*log(y0_d.^2 + S^2/4 + x0_d.*S + x0_d.^2));
u_d = -(gamma/(2*pi))*(atan((S/2 - x0_d)./y0_d) - ...
    atan((-S/2 - x0_d)./y0_d));

% Rotate coodinates back
u   = u_d*cos(phi) - v_d*sin(phi);
v   = v_d*cos(phi) + u_d*sin(phi);