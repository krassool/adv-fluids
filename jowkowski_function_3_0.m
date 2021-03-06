function sorted_panels = jowkowski_function_3_0(alpha_deg)

%% Set initial conditions
a   =  1       ; % whats this?
c   =  0.95    ; % Jowcow parameter
x_s = -0.04875 ; % x-equivilant thing
y_s =  0.05*1i  ; % y-equivilant thing

%% overwrite ics
a=1;
c=0.00;%0.95;
x_s=0.00%-0.04875;
y_s=0.0;%0.05*i;
alpha_deg=0; %degrees

%% Create physical variables

alpha = alpha_deg / 180*pi     ; % Angle of attack in radians
th    = 0 : pi/50 : 2*pi       ; % Theta for circle definition
z     = a*(cos(th)+1i*sin(th))  ; % Define the cylinder in the flow
z_c   = c*(cos(th)+1i*sin(th))  ; % Define the smaller transform circle (only for visuals)

%% Transformation 0 and 1
% We shift the circle with radius a relative to a Jowkowski 
% circle of radius c by the amounts xs and ys, such that the two circles 
% intersect on the x axis (to obtain the sharp trailing edge cusp).

w0 = z + x_s +y_s  ; 
w1 = w0+c^2./w0    ; % Apply the Joukowski transformation

%% Transformation 2, add the angle of attack ?

w2=w1*exp(-alpha*1i);

airfoil_panels=[real(w2);imag(w2)].';


n=length(airfoil_panels) ; % Number of panels to interate through
all_panels = zeros(n,4)  ; % Initialise the panels

for j=1:n
    if j==n % at the end, replace with the first ones
        all_panels(j,:) = [airfoil_panels(j,:),airfoil_panels(1,:)];
    else
        all_panels(j,:) = [airfoil_panels(j,:),airfoil_panels(j+1,:)];
    end
end

flip_panels   = (all_panels)          ; % Flip the top with the bottom%flip()
sorted_panels = circshift(flip_panels,-2) ; % Shift everything around by 2
% .. not sure why its 2?

end