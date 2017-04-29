% MCEN90018: Advanced Fluid Dynamics - Assignment 2
% ------------------------------------------------------------------------
% Mischka Kamener  539030                           Last modified: 28/4/16
%
% Calculates the circulation (gamma) of each of the vortex panels that make
% up the airfoil defined in x and y, in a uniform flow [U_inf, V_inf], such
% that the midpoints of each panel lie along a streamline. Additionaly, the
% first panel and last panel are set to have the same circulation, so that
% the Kutta condition is met.
function gamma = get_vortex_strengths(x, y, U_inf, V_inf)

% Get the number of panels from the coordinate matrices.
n_panels = length(x) - 1;

% Initialize matrices and vectors.
I = zeros(n_panels + 1);
V = zeros(n_panels + 1, 1);

% Set final row of I to satisfy Kutta condition, and final column to 
% satisfy constant streamfunction condition.
I(end,1) = 1;           % Kutta
I(end, end-1) = 1;      % Kutta
I(1:(end-1),end) = 1;   % Constant streamfunction

% Loop through all panels
for i = 1:n_panels
    
    % Set free stream streamfunction influence
    xm = mean(x(i:(i+1)));
    ym = mean(y(i:(i+1)));
    V(i) =  -U_inf*ym + V_inf*xm;
    
    % Find influence from all other panels
    for j = 1:n_panels
        
        % Set influence of panel j on midpoint of panel i;
        I(i,j) = induced_panel_streamfunction( ...
            x(i:(i+1)), y(i:(i+1)), x(j:(j+1)), y(j:(j+1)));
    end
end
% Solve for source panel strengths. Discard the streamfunction constant.
gamma = I\V;
gamma = gamma(1:(end-1));