%% Question 1: Source panel cylinder crossflow
% -- Define constants
n_panels = 8;
r = 1;
U_inf = 1;

% -- Create the panels and grid
% Streamlines go outside of [-2 2] range and return, so create larger grid 
[xp, yp] = meshgrid(-3:.01:3, -3:0.01:3);
panel_sep = 2*pi/n_panels;
theta    = flip((-pi+0.5*panel_sep):panel_sep:(pi+0.5*panel_sep));
x        = r*cos(theta);
y        = r*sin(theta);

% -- Solve for panel strengths and velocity field
q = get_source_strengths(x, y, U_inf);

% Initialsie velocity field
u = U_inf*ones(size(xp));
v = zeros(size(xp));

% Add velocity field contribution from each panel
for i = 1:n_panels
    [ui, vi] = source_panel_field(q(i), x(i:(i+1)), y(i:(i+1)), xp, yp);
    u = u + ui;
    v = v + vi;
end

% -- Plot field
% Find points withing cylinder and set velocity to NaN as area is not
% valid. (Streamlines won't pass through cylinder and contourf won't scale
% based on internal region)
[in, on] = inpolygon(xp,yp,x,y);
u(in & ~on) = NaN;
v(in & ~on) = NaN;

% Set the start of streamlines
starty = -2:0.25:2;
startx = -3*ones(size(starty));

% Create plot
figure
hold on
contourf(xp, yp, u, 20, 'LineStyle', 'none');
colormap winter
h = streamline(xp, yp, u, v, startx, starty);
set(h, 'Color', 'k');
fill(x, y, 'k');
xlabel('x')
ylabel('y')
title('8 Source Panel Cylinder')
axis equal
axis([-3 3 -2 2])
set(gcf, 'Position', [86 194 762 441]);

%% Question 2: Potential flow pressure vs panel method
% -- Find pressure at mid-points of panel
Uc = zeros(n_panels, 1);
for i = 1:n_panels
    
    % Find U_inf contribution
    phi = atan2(diff(y(i:(i+1))), diff(x(i:(i+1))));
    Uc(i) = U_inf*cos(2*pi - phi);
    
    % Find contributions from other panels
    for j = 1:n_panels
        % Tangential self influence is zero for source panel
        if (i == j), continue; end
        
        % Set self influence of panel j on panel i;
        [u_i, ~] = induced_panel_velocity(x(i:(i+1)), y(i:(i+1))...
            , x(j:(j+1)), y(j:(j+1)));
        Uc(i) = Uc(i) + q(j)*u_i;
    end
end
% Calculate pressure coefficient
Cp = 1-(Uc./U_inf).^2;

% Calculate theta at midpoints of panels
thetam = zeros(1, n_panels);
for i = 1:n_panels
    thetam(i) = mean(theta(i:(i+1)));
end

% -- Calculate theoretical pressure
a = 1;
theta_th = -pi:0.01:pi;

% Calculate circle coordinates.
x_th = a*cos(theta_th);  y_th = a*sin(theta_th);

% Calculate velocities from potential flow calcualtions.
u_th = U_inf*(1-a*(x_th.^2-y_th.^2)./(x_th.^2 + y_th.^2).^2);
v_th = -(2*a*U_inf.*x_th.*y_th)./(x_th.^2 + y_th.^2).^2;

% Calculate pressure coefficient.
Cp_th = 1-(sqrt(u_th.^2 + v_th.^2)/U_inf).^2;

% Plot pressures
figure
hold on
plot(thetam, Cp, 'rx');
plot(theta_th, Cp_th, 'b');
xlabel('\theta (rad)');
ylabel('C_{p}');
title('Panel Method Pressure vs Potential Flow Pressure')
legend({'Panel Midpoint', 'Theoretical'});
axis([-pi pi -3.1 1.1])
set(gcf, 'Position', [86 194 762 441]);

%% Question 4: Source panels with Jowkowski airfoil
% -- Create airfoil using Jowkowski transformation
aoa = 10*(pi/180); % Angle of attack (rad)
a = 1;          c = 0.95;
x_s = -0.0498;  y_s = 0.02;
n_panels = 200;
% Streamlines go outside of [-3 3] range and return, so create larger grid. 
[xp, yp] = meshgrid(-5:0.01:5, -3.5:0.01:3.5);
U_inf = 1;

% Determine angle shift so that panels start at trailing edge.
shift = atan(y_s./(c-x_s));
theta = flip((0-shift):2*pi/n_panels:(2*pi-shift));

% Transform and rotate coordinates to get airfoil.
z_cs = (a*cos(theta) + x_s) + 1i*(a*sin(theta) + y_s);
z_j  = exp(-1i*aoa)*(z_cs + (c^2./z_cs));
x = real(z_j);
y = imag(z_j);

% -- Use source panel method to get flow around airfoil.
% Solve for panel strengths
q = get_source_strengths(x, y, U_inf);

% Initialsie velocity field
u = U_inf*ones(size(xp));
v = zeros(size(xp));

% Add velocity field contribution from each panel
for i = 1:n_panels
    [ui, vi] = source_panel_field(q(i), x(i:(i+1)), y(i:(i+1)), xp, yp);
    u = u + ui;
    v = v + vi;
end

% Discard velocities within airfoil.
[in, on] = inpolygon(xp,yp,x,y);
u(in & ~on) = NaN;
v(in & ~on) = NaN;

% Set start of streamlines, streamline from trailing edge starts at x(1)
% y(1) and is added individually in the plot.
starty = -3:0.3:3;
startx = -5*ones(size(starty));

% Plot result
figure
hold on
contourf(xp, yp, u, 100, 'LineStyle', 'none');
caxis([0.5 1.5])
colormap(jet)
h1 = streamline(xp, yp, u, v, startx, starty);
h2 = streamline(xp, yp, u, v, x(1), y(1));
set(h1, 'Color', 'k');
set(h2, 'Color', 'r');
fill(x, y, 'k');
xlabel('x')
ylabel('y')
title('Source Panel Airfoil')
axis equal
axis([-5 5 -3 3])
set(gcf, 'Position', [86 194 762 441]);

%% Question 5: Vortex panels with Jowkowski airfoil
% Create airfoil using Jowkowski transformation
aoa = 10*(pi/180); % Angle of attack (rad)
a = 1;          c = 0.95;
x_s = -0.04875;  y_s = 0.05;
n_panels = 64;
% Streamlines go outside of [-3 3] range and return, so create larger grid. 
[xp, yp] = meshgrid(-5:0.01:5.2, -3.5:0.01:3.5);
U_inf = 1;

% Determine angle shift so that panels start at trailing edge.
shift = atan(y_s./(c-x_s));
theta = flip((0-shift):2*pi/n_panels:(2*pi-shift));

% Transform and rotate coordinates to get airfoil.
z_cs = (a*cos(theta) + x_s) + 1i*(a*sin(theta) + y_s);
z_j  = exp(-1i*aoa)*(z_cs + (c^2./z_cs));
x = real(z_j);
y = imag(z_j);

% -- Use vortex panel method to get flow around airfoil.
% Solve for panel circulations. Free stream flow only has U_inf component,
% V_inf = 0.
gamma = get_vortex_strengths(x, y, U_inf, 0);

% Initialsie velocity field
u = U_inf*ones(size(xp));
v = zeros(size(xp));

% Add velocity field contribution from each panel
for i = 1:n_panels
    [ui, vi] = vortex_panel_field( ...
        gamma(i), x(i:(i+1)), y(i:(i+1)), xp, yp);
    u = u + ui;
    v = v + vi;
end

% Discard velocities within airfoil.
[in, on] = inpolygon(xp,yp,x,y);
u(in & ~on) = NaN;
v(in & ~on) = NaN;

% Set start of streamlines, streamline from trailing edge starts at x(1)
% y(1) and is added individually in the plot.
starty = -3:0.3:3;
startx = -5*ones(size(starty));

figure
hold on
contourf(xp, yp, u, 100, 'LineStyle', 'none');
caxis([0.5 1.5])
colormap(jet)
xy_panel = stream2(xp, yp, u, v, [startx x(1)], [starty y(1)]);
h1 = streamline(xp, yp, u, v, startx, starty);
h2 = streamline(xp, yp, u, v, x(1), y(1));
set(h1, 'Color', 'k');
set(h2, 'Color', 'r');
fill(x, y, 'k');
xlabel('x')
ylabel('y')
title('Vortex Panel Airfoil')
axis equal
axis([-5 5 -3 3])
set(gcf, 'Position', [86 194 762 441]);

%% Question 7: Plot potential flow streamlines
% Determine starting points in untransformed plane
aoa = 10*(pi/180); % Angle of attack (rad)
a = 1;          c = 0.95;
x_s = -0.0498;  y_s = 0.02;

% Find potential flow solution to flow over cylinder.
U_inf = 1;
Gamma = 4*pi*U_inf*a*sin(asin(y_s/a)+aoa);
[x,y] = meshgrid(-6:0.01:12,-6:0.01:6);
z = x +1i.*y;
% Complex potential function for flow over cylinder is given by:
% w = U_inf*(z+a.^2./z)+(1i*Gamma./(2*pi)).*log(z./a);
dwdz = -a.^2*U_inf./(z.^2)+U_inf+1i*Gamma./(2*pi.*z);
u =  real(dwdz);
v = -imag(dwdz);
figure;
pcolor(x, y, real(sqrt(u.^2+v.^2))) ; shading flat ; colormap jet

%%

% -- Find streamlines for flow over cylinder and use Jowkowski transform to
%    get back to flow over airfoil.
% Set start locations based of Q5 specification, and reverse Jowkowski
% transform.
starty = -3:0.3:3;    startx = -5*ones(size(starty));
startz = startx + 1i*starty;
startz = [startz (1.8712 - 1i*0.32996)];

startz = startz*exp(1i*aoa);
 % Choose branch that gives solutions outside of cylider.
startz = 0.5*(startz-sqrt(startz.^2-4*c));
startz = startz - (x_s + 1i*y_s);
startz = startz*exp(-1i*aoa);

% Find cylinder boundary, and discard velocity values that lie within.
theta = 0:0.01:2*pi;
x_circ = a*cos(theta);
y_circ = a*sin(theta);
z_circ = x_circ + 1i*y_circ;
[in, on] = inpolygon(x,y,x_circ,y_circ);
u(in & ~on) = NaN;
v(in & ~on) = NaN;

% Find streamlines around cylinder, and convert to complex coordinates.
xy = stream2(x, y, u, v, real(startz), imag(startz));
z_pflow = zeros(10000, size(xy,2));
for i = 1:size(xy,2)
    xy_i = cell2mat(xy(i));
    z_pflow(:,i) = xy_i(:,1) + 1i*xy_i(:,2);
end

% Jowkowski transform the coordinates of the streamlines to get streamlines
% over airfoil.
z_pflow = z_pflow*exp(1i*aoa);
z_pflow = z_pflow + (x_s + 1i*y_s);
z_pflow = z_pflow + c^2./(z_pflow);
z_pflow = z_pflow*exp(-1i*aoa);

% Transform circle to plot airfoil.
z_circ = z_circ*exp(1i*aoa);
z_circ = z_circ + (x_s + 1i*y_s);
z_circ = z_circ + c^2./(z_circ);
z_circ = z_circ*exp(-1i*aoa);

figure
hold on
h1 = streamline(xy_panel);
set(h1, 'Color', 'k');
h2 = plot(real(z_pflow), imag(z_pflow), '--c');
fill(real(z_circ),imag(z_circ),'k');
xlabel('x')
ylabel('y')
title('Panel Method Streamlines vs Jowkowski Transform of Cylinder')
legend([h1(1) h2(1)], {'Panel Method', 'Jowkowski'})
axis equal
axis([-5 5 -3 3])
set(gcf, 'Position', [86 194 762 441]);