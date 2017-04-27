% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 8 Panel Cylinder
%
% Calculates the influence coefficients given thier midpoints
% Modified from the snippet code written by N.Hutchins for MCEN90018 

function I = panel_source_strength_1_0( Xi, Yi, Xj, Yj)

Xmi   = 0.5*(Xi(2)+Xi(1)); % midpoint of panel i
Ymi   = 0.5*(Yi(2)+Yi(1));
Phi_i = atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); % phi_i (eqn 24)

Xmj   = 0.5*(Xj(2)+Xj(1)); % midpoint of panel j
Ymj   = 0.5*(Yj(2)+Yj(1));

Phi_j = atan2((Yj(2)-Yj(1)),(Xj(2)-Xj(1))); %phi_j (eqn 23)

rij   = sqrt((Xmj - Xmi).^2 + (Ymj - Ymi).^2); % (eqn 22)

beta  = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
omega = beta - Phi_j; % (eqn 26)

x0p = rij.*cos(omega); % (eqn 27)
y0p = rij.*sin(omega); % (eqn 28)

S = sqrt((Xj(2) - Xj(1)).^2 + (Yj(2) - Yj(1)).^2);
a =-S/2; b =S/2; %(eqn 11)

vj = (1./(2*pi)).*(atan(((S./2)-x0p)./y0p)...
    -atan((-(S./2) - x0p)./y0p)); % eqn(30)

uj = (1./(2*pi)).*((-log((y0p.^2+((S.^2)./4)- (S.*x0p)+x0p.^2))./2)...
    + (log((y0p.^2 + ((S.^2)./4) + (S.*x0p) + x0p.^2))./2)); % eqn(29)

vi = uj.*sin(Phi_j-Phi_i)+vj.*cos(Phi_j-Phi_i); % eqn(31)

I = vi; % we do this as we set q2 to unity for now.