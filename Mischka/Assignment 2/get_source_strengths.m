% MCEN90018: Advanced Fluid Dynamics - Assignment 2
% ------------------------------------------------------------------------
% Mischka Kamener  539030                           Last modified: 28/4/16
%
% Calculates the source strength q of each of the panels that make up the 
% polygon defined in x and y, in a uniform flow U_inf, such that there is 
% no flow normal to the midpoint of each panel.
function q = get_source_strengths(x, y, U_inf)

% Get the number of panels from the coordinate matrices
n_panels = length(x) - 1;

% Initialize matrices and vectors
I = 0.5*ones(n_panels); % Initialize to self influence = (1/2)*q
V = zeros(n_panels, 1);

% Loop through all panels
for i = 1:n_panels
    
    % Set free stream velocity influence
    phi  = atan2(diff(y(i:(i+1))), diff(x(i:(i+1))));
    V(i) =  -U_inf*sin(2*pi - phi);
    
    % Find influence from all other panels
    for j = 1:n_panels
        
        % Self influence = (1/2)*q has already been set
        if (i == j), continue; end
        
        % Find the velocity influence normal to panel i due to panel j.
        [~, I(i,j)] = induced_panel_velocity(x(i:(i+1)), y(i:(i+1))...
            , x(j:(j+1)), y(j:(j+1)));
    end
end

% Solve for source panel strengths
q = I\V;