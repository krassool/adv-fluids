% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 'n' Panel Cylinder
%
% Estimates the flow field around an 'n' panel cylinder
%% Clear MATLAB environment, set format
clc , clear , close all %, format bank 

%% Create the panels and find the influsence co-efficients 

% Create airfoil panels using jowkowski 
aoa_degrees = 0 ;                                   % Angle of attack in degrees
panels      = jowkowski_function_2_0(aoa_degrees) ; % Create an airfoil in panels 
n_pan       = length(panels);                       % Number of panels
I = (zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ;   % Initialise influence 

% Calculate Influence
for m = 1 : n_pan ; % Loop throught each panel
    Xi = [panels(m,1),panels(m,3)]; % endpoints of panel j in x and y
    Yi = [panels(m,2),panels(m,4)];
    
    Phi_i(m) = -atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); % phi_i (eqn 24)
    
    for k=1:n_pan ; % Calculate the influence coeff on every other panel    
        Xj=[panels(k,1),panels(k,3)]; % Midpoints of panel i in x and y
        Yj=[panels(k,2),panels(k,4)];
        
        I(m,k)=panel_source_strength_1_0(Xi, Yi, Xj, Yj); % Find coeff
    end
end

% I=real(I);
I(eye(size(I))~=0) = 0.5;  % Where i==j hard code 0.5 strength (using logicals)
U_inf   = 1                     ;
V_inf_i = -U_inf*sin(2*pi-Phi_i); % find V_inf, flowing from left to right

q = I\V_inf_i; % Solve for source strength densities (q)

%% Find veloctities
tic
mesh_res = 0.01 ; % Meshgrid density (resolution for results)
[xp, yp] = meshgrid( -3:mesh_res:3 , -2:mesh_res:2 );
[u_hat,v_hat] = deal(zeros(size(xp))) ; % Initialise cartesian velocity directions 

% This next loop runs through each of the panels and sums the velocity
% contribution at each point in space as a result of the panels.

for n=1:n_pan ; % for each panel
    
    Xj=[panels(n,1),panels(n,3)];
    Yj=[panels(n,2),panels(n,4)];
    
    [u,v] = source_panel_on_point_vel( Xj , Yj , q(n) , xp , yp );
    
    u_hat=u_hat + u;
    v_hat=v_hat + v;
end
u_hat_inf = u_hat + U_inf;
time_pattern = toc

%% Caclulate the streamlines

% Set up simulation conditions
t0   = 0.00  ; % Initial time
tf   = 8.00  ; % Final time
h    = 0.01  ; % Step size
y_lim= 0.5;


y_range = (-y_lim:.2:y_lim).'   ;
% y_range = (-0.2:.1:0.2).';
ic0  = [ -3*ones(length(y_range),1) , y_range ];

% Initial conditions for streamlines
xs = ic0(:,1) ;
ys = ic0(:,2) ;

% Calculate streamlines in same fashion as fluids 1 
tic ; [xr, yr] = approx_streamline2(xs, ys, tf-t0, h, @flow_general , q , panels,U_inf);
time_streams   = toc

%% Plot results and make pretty
close all 
figure ; hold on ;

Xi=[panels(:,1),panels(:,3)]; % endpoints of panel j in x and y
Yi=[panels(:,2),panels(:,4)];

% Plot approximated cylinder with velocity field
plot(Xi, Yi, 'b-', 'LineWidth', 2.5) ; % Plot approximated shape
pcolor(xp, yp, real(sqrt(u_hat_inf.^2+v_hat.^2))) ; shading flat ; colormap jet
fill(panels(:,1),panels(:,2),[255 105 180]./256) ; % HOT PINK cylinder

% Create stream-lines
plot(xr.', yr.', 'r') ;

% Label plot and add features accordingly
axis equal  ; axis([-2 2 -2 2]) ; units = colorbar ;
xlabel(units,'m per s') ; xlabel('x (m)') ; ylabel('y (m)') ; 
legend('Streamlines')  ;
title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ') ;

%% Expose error in solution. Definition of airfoil known to be wrong

% Look at this, does make sense that its slow about the centre of the foil
v_mag = real(sqrt(u_hat_inf.^2+v_hat.^2)) ; 
figure ; spy(v_mag>=.5) ; title('Spying elements with velocity magnitude >= 0.5 m/s')
figure ; spy(v_mag<=.1) ; title('Spying elements with velocity magnitude <= 0.1 m/s')

exes = [panels(1:2,1) , panels(1:2,3) ] ;
whys = [panels(1:2,2) , panels(1:2,4) ] ;

figure ; plot(exes , whys) ; legend('line1','line2')
title('Showing first two panels')

% %% Plot results and make pretty
% 
% hold on ;  plot(xr.', yr.', 'k');
% 
% % quiver(xr(:,1), yr(:,1), xr(:,2)-xr(:,1), yr(:,2)-yr(:,1));
% axis equal  ; %axis([-2 2 -2 2])
% axis([-3 3 -3 3])
% xlabel('x (m)') ; ylabel('y (m)') ; legend('Streamlines')   ;
% title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ');
