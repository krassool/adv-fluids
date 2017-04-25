%jowkowski
close all

%% set initial conditions
a=1;
c=0.95;
x_s=-0.04875;
y_s=0.05*i;
alpha_deg=20; %degrees


%% create physical variables
alpha=alpha_deg/180*pi %angle of attack in radians
th = 0:pi/50:2*pi;      %theta for circle definition

%define the cylinder in the flow
z = a*(cos(th)+i*sin(th));

% Define the smaller transform circle (only for visuals)
z_c = c*(cos(th)+i*sin(th));

%% Transformation 1 
%We shift the circle with radius a relative to a Jowkowski 
%circle of radius c by the amounts xs and ys, such that the two circles 
%intersect on the x axis (to obtain the sharp trailing edge cusp).

w0 =z + x_s +y_s;

figure
hold on
plot(real(z),imag(z),'b')
plot(real(w0),imag(w0),'g')
axis([-1.1 1.1 -1.1 1.1])
daspect([1 1 1])


%% Transformation 2 
%We apply the Joukowski transformation

w1=w0+c^2./w0;

%% Transformation 2 
%We add the angle of attack ?

w2=w1*exp(-alpha*i);


figure
plot(real(w2),imag(w2),'g.')
axis([-5 5 -3 3])
daspect([1 1 1])

airfoil_panels=[real(w2),imag(w2)];
