% Test Vortex Panel Script:
clc, clear, close all

n_el    = 2 
panel_x = linspace(-4,4,n_el+1)
panel_y = zeros(1,length(panel_x))

n = length(panel_x)

% Create end-points
for i=1:n
    if i==n % at the end, replace with the first ones
        all_panels(i,:) = [new_panels(i,:),new_panels(1,:)];
    else
        all_panels(i,:) = [new_panels(i,:),new_panels(i+1,:)];
    end
end


    
% Calculate influence
for m = 1 : n_pan  % Loop throught each panel
    Xi  = [panels(m,1),panels(m,3)] ; % end?points of panel j in x and y
    Yi  = [panels(m,2),panels(m,4)] ;
    
    Phi_i(m) = atan2((Yi(2) -Yi(1)),(Xi(2) - Xi(1))); % Phi_i (eqn 24) 
    
    for k=1 : n_pan   % Calculate the influence coeff on every other panel    
        Xj=[panels(k,1),panels(k,3)]; % Midpoints of panel i in x and y
        Yj=[panels(k,2),panels(k,4)];
        
        I(m,k) = vortex_panel_strength_2_0(Xi, Yi, Xj, Yj); % Find coeff
    end
end

I(eye(size(I))~=0) = 0.5;  % Where i==j hard code 0.5 strength (using logicals)
I_append=zeros(1,size(I,1));
I_append(1)=1; I_append(end)=1;

I_con=[I;I_append];

V_inf_i     =  -U_inf*sin(2*pi-Phi_i) % find V_inf, flowing from left to right
V_inf_i_con = [ V_inf_i ; 1] ;
% gam         = I_con\V_inf_i_con       % Solve for source strength densities (q)

gam         = I\V_inf_i       % Solve for source strength densities (q)
% gam=gam/n_pan;
%% Find veloctities
tic
mesh_res      = 0.01 ; % Meshgrid density (resolution for results)
[xp, yp]      = meshgrid( -2:mesh_res:2 , -2:mesh_res:2 ) ; 
[u_hat,v_hat] = deal(zeros(size(xp))) ; % Initialise cartesian velocity directions 

% This next loop runs through each of the panels and sums the velocity
% contribution at each point in space as a result of the panels.

for n=1:n_pan ; % for each panel
    
    Xj = [panels(n,1),panels(n,3)] ;
    Yj = [panels(n,2),panels(n,4)] ;
    
    [u,v] = vor_panel_on_point_vel( Xj , Yj , gam(n) , xp , yp );
    
    u_hat=u_hat + u ;
    v_hat=v_hat + v ;
end

u_hat_inf = u_hat + U_inf;
time_pattern = toc