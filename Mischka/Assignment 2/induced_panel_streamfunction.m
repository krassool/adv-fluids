% MCEN90018: Advanced Fluid Dynamics - Assignment 2
% ------------------------------------------------------------------------
% Mischka Kamener  539030                           Last modified: 28/4/16
%
% Function based of "snippet.m" written by Nick Hutchins. Calculates the
% streamfunction by the vortex panel j on the midpoint of panel i.
function psi_i = induced_panel_streamfunction(xi, yi, xj, yj)

% Find midpoints of panel i and j.
xmi   = mean(xi);
ymi   = mean(yi);
xmj   = mean(xj);
ymj   = mean(yj);

% Calculate distance between mipoints.
r = sqrt((xmj - xmi).^2 + (ymj - ymi).^2);

% Calculate angles.
phi_j =atan2(diff(yj),diff(xj));
beta = atan2((ymi - ymj),(xmi - xmj));
omega = beta - phi_j;

% Calculate x, y distance between panel midpoints in frame of panel j.
x0_d = r.*cos(omega);
y0_d = r.*sin(omega);
S = sqrt(diff(xj).^2 + diff(yj).^2); % Length of panel j.
a = -S/2;
b = S/2;

% Calculate streamfunction psi at the midpoint of panel i. Do not need to
% change coordinates, as streamfunction is a scalar. gamma_j = 1, as we 
% take it out of this equation. 
psi_i = (1./(2.*pi)).*((((x0_d-b)./2).*log((x0_d-b).^2 + y0_d.^2) + ...
    y0_d.*atan((x0_d-b)./y0_d) + b) - ...
    (((x0_d-a)./2).*log((x0_d-a).^2+y0_d.^2) + ...
    y0_d.*atan((x0_d-a)./y0_d) + a));


end

