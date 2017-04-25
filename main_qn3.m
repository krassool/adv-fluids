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
% panels = create8pan
I=(zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ; % Initialise influence 

% Calculate influence
figure ; hold on ;

for m=1:n_pan; % Loop throught each panel
    Xi=[panels(m,1),panels(m,3)]; % end?points of panel j in x and y
    Yi=[panels(m,2),panels(m,4)];
    
    Phi_i(m)=atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); % phi_i (eqn 24)
    plot(Xi, Yi, 'b-', 'LineWidth', 2.5) ;          % Plot approximated cylinder
    
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
    
    [u,v] = flow_field_cyl_1_0( Xj , Yj , q(n) , xp , yp );
    
    u_hat=u_hat + u;
    v_hat=v_hat + v;
    
    [u_surf,v_surf] = flow_field_cyl_1_0( Xj , Yj , q(n) , Xmj , Ymj );
    
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

Cp=1-(U_c/U_inf).^2;

quiver(Xmj,Ymj,U_c_u,U_c_v);



%% Plot panel estimation

h = pcolor(xp, yp, sqrt(u_hat_inf.^2+v_hat.^2)) ;
set(h, 'EdgeColor', 'none') ; colormap jet
% quiver(xp, yp, u_hat_inf, v_hat, 2, 'k')
fill(panels(:,1),panels(:,2),[255 105 180]./256)

%% Plot the streamlines

% Set up simulation conditions
t0   = 0     ; % Initial time
tf   = 6.00  ; % Final time
h    = 0.01  ; % Step size

y_range = (-2:.25:2).';
ic0  = [ -3*ones(length(y_range),1) , y_range ];

% Initial conditions for streamlines
xs = ic0(:,1) ;
ys = ic0(:,2) ;
c = [n,q.']   ; % Pack up concstants matrix 

% Calculate streamlines in same fashion as fluids 1 
tic ; [xr, yr] = approx_streamline2(xs, ys, tf-t0, h, @flow_general , q , panels, U_inf);
time_streams = toc
%% Plotting fucntions
% Plot results and make pretty
hold on ;  plot(xr.', yr.', 'b') ;
quiver(xr(:,1), yr(:,1), xr(:,2)-xr(:,1), yr(:,2)-yr(:,1));
axis equal  ; axis([-2 2 -2 2]) ;  colorbar ; 
% axis([-4 4 -4 4])
xlabel('x (m)') ; ylabel('y (m)') ; legend('Streamlines')    ;

title('Flow over and 8 Panel Cylinder (w.page, k.rassool) ') ;