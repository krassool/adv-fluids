% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 8 Panel Cylinder
%
% Streamline solver using RK4 and example_code_figure5 from N.Hutchins

function [state] = flow_general(state0, h, n, q , panels, u_inf)

% Create space for new solution
state = [state0,zeros(length(state0),n)] ;

% Iterate over number of steps
for i = 1:n
    
    % Determine RK4 parameters
    k1 = get_velocities(state(:,i), q , panels, u_inf);
    k2 = get_velocities(state(:,i) + (h/2)*k1, q , panels, u_inf);
    k3 = get_velocities(state(:,i) + (h/2)*k2, q , panels, u_inf);
    k4 = get_velocities(state(:,i) +  h*k3   , q , panels, u_inf);
    
    % Update RK4 estimate for next step
    state(:,i+1) = state(:,i) + (h/6)*(k1+2*k2+2*k3+k4);
end

function state_derivative = get_velocities( state , q , panels , u_inf)

x = state(1) ; y = state(2); % Previous state definition

% Get velocity contributions from each panel
[u,v] = deal(0);

for n  = 1:length(panels)
    Xj = [panels(n,1),panels(n,3)] ; % Panel endpoints in X and Y
    Yj = [panels(n,2),panels(n,4)] ;
    
    [u_tmp,v_tmp] = flow_field_cyl_1_0( Xj , Yj , q(n) , x , y );    
    u = u+u_tmp;
    v = v+v_tmp;
end
u = u+u_inf ;

state_derivative = [u ; v]; 