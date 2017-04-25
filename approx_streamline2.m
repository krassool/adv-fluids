% William Page (587000) - Kevin Rassool (xxxxxx) ;
% Semester 2 2015 - University of Melbourne        ; Started:     21/4/17
% MCEN90018 - Advanced Fluid Dynamics              ; Last Edited: 29/4/17
% Assignment 2 : Panel Methods - 8 Panel Cylinder
%
% Plots the streamlines around an 8 panel cylinder

function [xr, yr] = approx_streamline2(xs, ys, tf, h, fn, c)

n = tf/h;    % Number of steps
[xr,yr] = deal([xs , zeros(length(xs),n)]) ; % Pre-allocate xr and yr array 

for i = 1:length(xs)
    state0 = [xs(i) ; ys(i)]; % Define state
    state = fn(state0, h, n, c); % Run simulation

    % Copy into return array
    xr(i,:) = state(1,:);   % x-coords in first row
    yr(i,:) = state(2,:);   % y-coords in second row
end
end