% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 'n' Panel Cylinder
%
% Estimates the flow field around an 'n' panel cylinder
%% Clear MATLAB environment, set format
clc , clear , close all %, format bank 

%% Create the panels and find the influsence co-efficients 
U_inf = 1 ;
n_pan = 8 ; % Number of panels to use
panels = n_panel_circle(n_pan) ;  % Define the panels for 8 (make mathematical l8r)

%streamline solution constants
k_val=pi*U_inf;

% panels = create8pan
I=(zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ; % Initialise influence 

% Calculate influence
figure ; hold on ;

for m=1:n_pan; % Loop throught each panel
    Xi=[panels(m,1),panels(m,3)]; % end?points of panel j in x and y
    Yi=[panels(m,2),panels(m,4)];
    
    Phi_i(m)=atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); % phi_i (eqn 24)
    
    for k=1:n_pan ; % Calculate the influence coeff on every other panel    
        Xj=[panels(k,1),panels(k,3)]; % Midpoints of panel i in x and y
        Yj=[panels(k,2),panels(k,4)];
        
        I(m,k)=panel_source_strength_1_0(Xi, Yi, Xj, Yj); % Find coeff
    end
end

I(eye(size(I))~=0) = 0.5;  % Where i==j hard code 0.5 strength (using logicals)

V_inf_i = -U_inf*sin(2*pi-Phi_i); % find V_inf, flowing from left to right

q = I\V_inf_i % Solve for source strength densities (q)

%% Find veloctities
tic

%midpoint matrix
r_mid=1.001;
panels_points = n_panel_circle_soft(n_pan,r_mid) 
Xmj=0.5*(panels_points(:,1)+panels_points(:,3))
Ymj=0.5*(panels_points(:,2)+panels_points(:,4)) 

mesh_res = 0.01 ; % Meshgrid density (resolution for results)
[xp, yp] = meshgrid( -2:mesh_res:2 , -2:mesh_res:2 );
[u_hat,v_hat] = deal(zeros(size(xp))) ; % Initialise cartesian velocity directions 
[u_hat_surf,v_hat_surf]=deal(zeros(size(Xmj))) ; % Initialise midpoint velocity directions 
% This next loop runs through each of the panels and sums the velocity
% contribution at each point in space as a result of the panels.   

for n=1:n_pan  % for each panel
    
    Xj=[panels(n,1),panels(n,3)];
    Yj=[panels(n,2),panels(n,4)];
    
    [u,v] = source_panel_on_point_vel( Xj , Yj , q(n) , xp , yp );
    
    u_hat=u_hat + u;
    v_hat=v_hat + v;
    
    [u_surf,v_surf] = source_panel_on_point_vel( Xj , Yj , q(n) , Xmj , Ymj );
    
    u_hat_surf=u_hat_surf + u_surf;
    v_hat_surf=v_hat_surf + v_surf;
end

u_hat_inf = u_hat + U_inf;
time_pattern = toc

%% Calculate Cp

%this is in XXX coordainate frame.
U_c_u=u_hat_surf + U_inf;
U_c_v=v_hat_surf;% + V_inf_i ;
U_c=sqrt(U_c_u.^2+U_c_v.^2);
%Cp opbtained from panel methods
Cp=1-(U_c/U_inf).^2;
%Cp obtained from Theory
th = 0:2*pi/n_pan:2*pi;
Cp_th=(1-4*sin(th).^2).';

figure
hold on
plot(Cp_th,'r*')
plot(Cp,'bx')


% Cp_matrix=[Cp Cp_th];
%% Plot the streamlines

% Set up simulation conditions
t0   = 0     ; % Initial time
tf   = 6  ; % Final time
h    = 0.01  ; % Step size

y_range = (-2:.1:2).';
ic0  = [ -3*ones(length(y_range),1) , y_range ];

% Initial conditions for streamlines
xs = ic0(:,1) ;
ys = ic0(:,2) ;

% Calculate streamlines in same fashion as fluids 1 
tic ; [xr, yr] = approx_streamline2(xs, ys, tf-t0, h, @flow_general , q , panels, U_inf);
[xr_s, yr_s] = approx_streamline2(xs, ys, tf-t0, h, @flow_general_streamline , k_val , panels, U_inf);
time_streams = toc
%% Plot results and make pretty
close all 
figure ; hold on ;

Xi=[panels(:,1),panels(:,3)]; % endpoints of panel j in x and y
Yi=[panels(:,2),panels(:,4)];

% Plot approximated cylinder with velocity field
plot(Xi, Yi, 'b-', 'LineWidth', 2.5) ; % Plot approximated cylinder
pcolor(xp, yp, sqrt(u_hat_inf.^2+v_hat.^2)) ; shading flat ; colormap jet
fill(panels(:,1),panels(:,2),[255 105 180]./256) ; % HOT PINK cylinder

% Create stream-lines
plot(xr.', yr.', 'b') ;

% Create theory stream-lines
plot(xr_s.', yr_s.', 'r') ;

% Label plot and add features accordingly
axis equal  ; axis([-2 2 -2 2]) ;  h = colorbar ;
xlabel(h,'m per s') ; xlabel('x (m)') ; ylabel('y (m)') ; 
legend('Streamlines')    ;
title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ') ;
% quivers(xr(:,100), yr(:,100), ((xr(:,101)-xr(:,100))/h), ((yr(:,101)-yr(:,100))/h) ... 
%     ,0.5,1,'m/s','k');