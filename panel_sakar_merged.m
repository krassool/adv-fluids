%  Panel Code in MATLAB
%
%  Coded by L. sankar, April 1997
%
clear
close all
%% Create the panels and find the influsence co-efficients
U_inf  = 1  ;
% n_pan  = 64 ; % Number of panels to use
% panels = n_panel_circle(n_pan) ;  % Define the number of approximation panels
% I      =(zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ; % Initialise influence

% Create airfoil panels using jowkowski
aoa_degrees = 00 ;                                   % Angle of attack in degrees
panels      = jowkowski_function_5_0(aoa_degrees) ; % Create an airfoil in panels
n_pan       = length(panels);                       % Number of panels
% n_pan       = 5;
% panels      = n_panel_circle(n_pan) ;  % Define the number of approximation panels
I           = (zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ;   % Initialise influence


% Read Angle of Attack
%
alpha = aoa_degrees/pi*180;
%
% read number of points on the upper side of airfoil
%
nu = n_pan/2
%
% read number of points on the lower side of airfoil
%
nl = n_pan/2
%
% read Flag that states if this airfoil is symmetric
% if isym > 0 then airfoil is assumed symmetric
%
isym = 1
%
% Read a scaling factor
% The airfoil y- ordinates will be multiplied by this factor
%
factor= 3;

if(isym>0)
  nl = nu;
end
%
% Allocate storage for x and y
%
x = zeros(1,100);
y = zeros(1,100);
%
% Read the points on the upper surface
%
numel(nu)
   x = panels(:,1);
   y = panels(:,2);
   
   x_append=x(1);
   y_append=y(1);
   
   x=[x;x_append];
   y=[y;y_append];
   
%
% If the airfoil is not symmetric, read lower side ordinates too..
%

%
% Plot the airfoil on window #1
%
plot(x,y);
n=ceil(length(x)/2);
A=zeros(n+1,n+1);
ds=zeros(1,n);
pi=4. * atan(1.0);
%
% Assemble the Influence Coefficient Matrix A
%
 for i = 1:n
   t1= x(i+1)-x(i);
   t2 = y(i+1)-y(i);
   ds(i) = sqrt(t1*t1+t2*t2);
 end
 
for j = 1:n
 a(j,n+1) = 1.0;
 for i = 1:n
   if i == j
	 
     fprintf('This is Ds : %f \n',ds(i))
     fprintf('This is the log bit of Ds : %f \n',(log(0.5*ds(i)) - 1.0) )
     
     a(i,i) = ds(i)/(2.*pi) *(log(0.5*ds(i)) - 1.0);
   else
     xm1 = 0.5 * (x(j)+x(j+1));
     ym1 = 0.5 * (y(j)+y(j+1));
     dx  = (x(i+1)-x(i))/ds(i);
     dy  = (y(i+1)-y(i))/ds(i);
     t1  = x(i) - xm1;
     t2  = y(i) - ym1;
     t3  = x(i+1) - xm1;
     t7  = y(i+1) - ym1;
     t4  = t1 * dx + t2 * dy;
     t5  = t3 * dx + t7 * dy;
     t6  = t2 * dx - t1 * dy;
     t1  = t5 * log(t5*t5+t6*t6) - t4 * log(t4*t4+t6*t6);
     t2  = atan2(t6,t4)-atan2(t6,t5);
     a(j,i) = (0.5 * t1-t5+t4+t6*t2)/(2.*pi);
   end
 end
a(n+1,1) = 1.0;
a(n+1,n) = 1.0;
end
%
% Assemble the Right hand Side of the Matrix system
%
rhs=zeros(n+1,1);
alpha = alpha * pi /180;
xmid=zeros(n,1);
for i = 1:n
  xmid(i,1) = 0.5 * (x(i) + x(i+1));
  ymid = 0.5 * (y(i) + y(i+1));
  rhs(i,1) = ymid * cos(alpha) - xmid(i) * sin(alpha);
end
gamma = zeros(n+1,1);
%
% Solve the syetm of equations
% In MATLAB this is easy!
%
gamma =   a\rhs;
gam = gamma;
I_matrix = a
V_inf_matrix=rhs;

%% Find veloctities
tic
mesh_res      = 0.01 ; % Meshgrid density (resolution for results)
[xp, yp]      = meshgrid( -2:mesh_res:2 , -2:mesh_res:2 ) ; 
[u_hat,v_hat] = deal(zeros(size(xp))) ; % Initialise cartesian velocity directions 

% This next loop runs through each of the panels and sums the velocity
% contribution at each point in space as a result of the panels.

for gg=1:n; % for each panel
    
    Xj = [panels(gg,1),panels(gg,3)] ;
    Yj = [panels(gg,2),panels(gg,4)] ;
    
    [u,v] = vor_panel_on_point_vel( Xj , Yj , gam(gg) , xp , yp );
    
    u_hat=u_hat + u ;
    v_hat=v_hat + v ;
end

u_hat_inf = u_hat + U_inf;
time_pattern = toc
%% Solve the streamlines

% Set up simulation conditions
t0   = 0     ; % Initial time
tf   = 6     ; % Final time
h    = 0.01  ; % Step size

y_range = (-2:.25:2).'; % Range over which to seen line for flow definition
ic0  = [ -3*ones(length(y_range),1) , y_range ]; % % Initial condition matrix

% Initial conditions for streamlines
xs = ic0(:,1) ;
ys = ic0(:,2) ;

% % Calculate streamlines in same fashion as fluids 
tic ; [xr, yr] = approx_streamline2(xs, ys, tf-t0, h, @flow_vor_general , gam , panels, U_inf);
time_streams = toc

%% Plot results and make pretty
close all 
figure ; hold on ;

Xi=[panels(:,1),panels(:,3)]; % endpoints of panel j in x and y
Yi=[panels(:,2),panels(:,4)];

% Plot approximated cylinder with velocity field
plot(Xi, Yi, 'b-', 'LineWidth', 2.5) ; % Plot approximated cylinder
pcolor(xp, yp, real(sqrt(u_hat_inf.^2+v_hat.^2))) ; shading flat ; colormap jet

% Create stream-lines
plot(xr.', yr.', 'r') ;
% Indicate streamline direction and magnitude
quivers(xr(:,100), yr(:,100), (xr(:,101)-xr(:,100))./h, (yr(:,101)-yr(:,100))./h , 0.5 , 1 , 'm/s' , 'k')

% Label plot and add features accordingly
fill(panels(:,1),panels(:,2),[255 105 180]./256) ; % HOT PINK cylinder
axis equal  ; axis([-2 2 -2 2]) ; units = colorbar ;
xlabel(units,'m per s') ; xlabel('x (m)') ; ylabel('y (m)') ; 
caxis([0 5]);
legend('Streamlines')    ;
title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ') ;