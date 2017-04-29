Xi = [0.3827    0.9239]; Yi = [0.9239    0.3827]; % coordinates of panel i
Xmi=0.5*(Xi(2)+Xi(1)); % midpoint of panel i
Ymi=0.5*(Yi(2)+Yi(1));
Phi_i=atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); %phi_i (eqn 24)

q_j = 1;

Xj = [-0.9239   -0.3827]; Yj = [0.3827    0.9239]; % coordinates of panel j
Xmj=0.5*(Xj(2)+Xj(1)); % midpoint of panel j
Ymj=0.5*(Yj(2)+Yj(1));
Phi_j=atan2((Yj(2)-Yj(1)),(Xj(2)-Xj(1))); %phi_j (eqn 23)

rij = sqrt((Xmj - Xmi).^2 + (Ymj - Ymi).^2); % (eqn 22)

beta = atan2((Ymi - Ymj),(Xmi - Xmj)); % (eqn 25)
omega = beta - Phi_j; % (eqn 26)

x0p = rij.*cos(omega); % (eqn 27)
y0p = rij.*sin(omega); % (eqn 28)

S = sqrt((Xj(2) - Xj(1)).^2 + (Yj(2) - Yj(1)).^2);
a =-S/2; b =S/2; %(eqn 11)

vj = (q_j./(2*pi)).*(atan(((S./2)-x0p)./y0p)...
    -atan((-(S./2) - x0p)./y0p)); % eqn(30)

uj = (q_j./(2*pi)).*((-log((y0p.^2+((S.^2)./4)- (S.*x0p)+x0p.^2))./2)...
    + (log((y0p.^2 + ((S.^2)./4) + (S.*x0p) + x0p.^2))./2)); % eqn(29)

vi=uj.*sin(Phi_j-Phi_i)+vj.*cos(Phi_j-Phi_i) % eqn(31)
