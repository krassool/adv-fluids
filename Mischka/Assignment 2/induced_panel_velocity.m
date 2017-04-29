% MCEN90018: Advanced Fluid Dynamics - Assignment 2
% ------------------------------------------------------------------------
% Mischka Kamener  539030                           Last modified: 28/4/16
%
% Function based of "snippet.m" written by Nick Hutchins. Calculates the
% induced velocities ui and vi in the frame of panel i, by the source 
% panel j on the midpoint of panel i.
function [ui, vi] = induced_panel_velocity(xi, yi, xj, yj)

% Find midpoints of panel i and j.
xmi   = mean(xi);
ymi   = mean(yi);
xmj   = mean(xj);
ymj   = mean(yj);

% Calculate distance between mipoints.
r = sqrt((xmj - xmi).^2 + (ymj - ymi).^2);

% Calculate angles.
phi_i =atan2(diff(yi), diff(xi));
phi_j =atan2(diff(yj),diff(xj));
beta = atan2((ymi - ymj),(xmi - xmj));
omega = beta - phi_j;

% Calculate x, y distance between panel midpoints in frame of panel j.
x0_d = r.*cos(omega);
y0_d = r.*sin(omega);
S = sqrt(diff(xj).^2 + diff(yj).^2); % Length of panel j

% Find velocities frame of panel j. q_j = 1, as we take it out of this
% equation.
vj = (1./(2*pi)).*(atan(((S./2)-x0_d)./y0_d)...
    -atan((-(S./2) - x0_d)./y0_d));
uj = (1./(2*pi)).*((-log((y0_d.^2+((S.^2)./4)- (S.*x0_d)+x0_d.^2))./2)...
    + (log((y0_d.^2 + ((S.^2)./4) + (S.*x0_d) + x0_d.^2))./2)); % eqn(29)

% Rotate velocities into frame of panel i.
ui = uj.*cos(phi_j-phi_i)-vj.*sin(phi_j-phi_i);
vi = uj.*sin(phi_j-phi_i)+vj.*cos(phi_j-phi_i);