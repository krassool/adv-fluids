%% Define constants
n_panels = 8;
r = 1;
U_inf  = 1;

%% Create the panels and grid
[xp, yp] = meshgrid(-3:.1:3, -2:0.1:2);
panel_sep = 2*pi/n_panels;
theta    = (-pi+0.5*panel_sep):panel_sep:(pi+0.5*panel_sep);
theta    = flip(theta);
x = r*cos(theta);
y = r*sin(theta);

plot(x,y)

q = get_source_strengths(x, y, U_inf)