% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 8 Panel Cylinder
%
% Streamline solver using RK4 and example_code_figure5 from N.Hutchins

function [state] = flow_general_streamline(state0, h, n, k_val, panels ,u_inf)
%state contains x and y values
%h is the timestep
%n is the number of steps to be taken

% Create space for new solution
state = [state0,zeros(length(state0),n)] ;

% Iterate over number of steps
for i = 1:n
    
    % Determine RK4 parameters
    k1 = get_velocities(state(:,i), k_val ,panels, u_inf);
    k2 = get_velocities(state(:,i) + (h/2)*k1, k_val ,panels, u_inf);
    k3 = get_velocities(state(:,i) + (h/2)*k2, k_val ,panels, u_inf);
    k4 = get_velocities(state(:,i) +  h*k3   , k_val ,panels, u_inf);
    
    % Update RK4 estimate for next step
    state(:,i+1) = state(:,i) + (h/6)*(k1+2*k2+2*k3+k4);
    
%     % Check for singularity
%     if sing_check(state(:,i), (h/6)*(k1+2*k2+2*k3+k4), q_all) == 1
%         state(:,i+1) = state(:,i);
%     end
end

function state_derivative = get_velocities( state , k , panels , u_inf)
x = state(1) ; y = state(2); % Previous state definition

u = u_inf+ (k*(x^2-y^2))/(pi*(x^2+y^2)^2) ;
v = -(2*k*y*x)/(pi*(x^2+y^2)^2);

state_derivative = [u ; v]; 

% function is_singularity = sing_check( state, state_derivative, c )
% 
% % Unpack constants
% ls = c(3); le = c(4);
% 
% % Define singularity locations
% singularity = [0 ls; 0 le; 0 -ls; 0 -le];
% 
% is_singularity = 0;
% % If the vector of the velocity is greater than distance to the singularity 
% for i = 1:length(singularity(:,1))
%     if norm(state.' - singularity(i,:)) < norm(state_derivative)
%         is_singularity = 1;
%     end
% end