% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 'n' Panel Cylinder
%
% Estimates the flow field around an 'n' panel cylinder
%% Clear MATLAB environment, set format
clear , %, format bank 

%% Physical Constants
aoa=10; %degrees
alpha   =aoa/180*pi; %angle of attack in radians
a   =  1       ;
c   =  0.95    ;
x_s =  -0.04875 ;
y_s =    0.05*1i ;
n_steps=64     ;
U_inf       = 1;

gam=0; , panels=0; %initialise crap;
%% Set up mesh
xtent   = 6;
mesh_res= 0.01 ; % Meshgrid density (resolution for results)
[x,y]   = meshgrid(-xtent:mesh_res:xtent,-6:0.01:6);

%% Find veloctities

z       = x + 1i.*y;
Tau = -4*pi*U_inf*a*sin(asin(y_s/a)+alpha);
dwdz_circul = (-1i*Tau)./(2*pi*(z) );
dwdz_cyl    = -a.^2*U_inf./(z.^2)+ U_inf;
dwdz=dwdz_cyl + dwdz_circul ;

dwdz = -a.^2*U_inf./(z.^2)+U_inf+1i*Tau./(2*pi.*z);


u =  real(dwdz);
v = -imag(dwdz);


%% Solve the streamlines

% Set up simulation conditions
t0   = 0     ; % Initial time
tf   = 16     ; % Final time
h    = 0.01  ; % Step size

y_range = (-2:.25:2).'; % Range over which to seen line for flow definition
ic0  = [ -xtent*ones(length(y_range),1) , y_range ]; % % Initial condition matrix

% Initial conditions for streamlines
xs = ic0(:,1) ;
ys = ic0(:,2) ;

% % % Calculate streamlines in same fashion as fluids 
tic ; [xr_1, yr_1] = approx_streamline2(xs, ys, tf-t0, h, @flow_vor_complex_circulation , gam , panels, U_inf);
time_streams = toc


%% convert streamlines to complex to conduct transofrm

z_1=xr_1+1i*yr_1;

z_2 = z_1*exp(1i*alpha);
z_3 = z_2 + (x_s + 1i*y_s);
z_4 = z_3 + c^2./(z_3);
z_5 = z_4*exp(-1i*alpha);

xr=real(z_5);
yr=imag(z_5);
%% Plot results and make pretty

hold on ;
% 
xr(7,6.4/h:end)=xr(7,6.4/h);
yr(7,6.4/h:end)=yr(7,6.4/h);

xr(6,7/h:end)=xr(6,7/h);
yr(6,7/h:end)=yr(6,7/h);

% Create stream-lines
hold on
plot(xr.', yr.', 'b') ;
panels      = jowkowski_function_5_0(aoa) ; % Create an airfoil in panels 
n_pan       = length(panels);                       % Number of panels
% panels      = n_panel_circle(n_pan) ;  % Define the number of approximation panels

 fill(panels(:,1),panels(:,2),[255 105 180]./256) ; % HOT PINK cylinder
axis equal  ; axis([-xtent xtent -2 2]) ; units = colorbar ;
xlabel(units,'m per s') ; xlabel('x (m)') ; ylabel('y (m)') ; 
legend('Streamlines')    ;  
title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ') ;

