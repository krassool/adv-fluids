%jowkowski
close all

%% set initial conditions
a=1;
c=0.95;
x_s=-0.04875;
y_s=0.05*i;
alpha_deg=0; %degrees


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

%% accounting for the panels

airfoil_panels=[real(w2);imag(w2)].';
n=length(airfoil_panels);
all_panels = zeros(n,4);

for j=1:n
    if j==n % at the end, replace with the first ones
        all_panels(j,:) = [airfoil_panels(j,:),airfoil_panels(1,:)];
    else
        all_panels(j,:) = [airfoil_panels(j,:),airfoil_panels(j+1,:)];
    end
end

flip_panels=flip(all_panels);
sorted_panels=circshift(flip_panels,-1);
clipped_panels=sorted_panels(1:end-1,:);


%% Plotting
figure;
   Xj = [clipped_panels(1,1),clipped_panels(1,3)]
   Yj = [clipped_panels(1,2),clipped_panels(1,4)]
   plot(Xj,Yj,'r*','LineWidth', 2.5)
 hold on ;
% for p=1:n
%    Xj = [clipped_panels(p,1),clipped_panels(p,3)]
%    Yj = [clipped_panels(p,2),clipped_panels(p,4)]
%    plot(Xj,Yj,'-b','LineWidth', 2.5)
%    pause
% end
