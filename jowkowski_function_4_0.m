function clipped_panels = jowkowski_function_4_0(alpha_deg)

%% set initial conditions
a   =  1       ;
c   =  0.95    ;
x_s =  0.04875 ;
y_s =  0.05*1i ;
n_steps=8     ;
%% set temp initial conditions
a=1;
c=0.5;%0.95;
x_s=0.0%-0.04875;
y_s=0.0;%0.05*i;
alpha_deg=0; %degrees

%% create physical variables
alpha   =alpha_deg/180*pi %angle of attack in radians
th      = linspace(0,2*pi,n_steps);      %theta for circle definition

%define the cylinder in the flow
z = a*(cos(th)+1i*sin(th));

% Define the smaller transform circle (only for visuals)
z_c = c*(cos(th)+1i*sin(th));

%% Transformation 1 
%We shift the circle with radius a relative to a Jowkowski 
%circle of radius c by the amounts xs and ys, such that the two circles 
%intersect on the x axis (to obtain the sharp trailing edge cusp).

w0 =z + x_s +y_s;

%% Transformation 2 
%We apply the Joukowski transformation

w1=w0+c^2./w0;

%% Transformation 2 
%We add the angle of attack ?

w2=w1*exp(-alpha*1i);

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

%% Plotting (optional)


%% Plotting

figure(1) ; hold on ;
for p=1:n-1
   Xj = [clipped_panels(p,1),clipped_panels(p,3)]
   Yj = [clipped_panels(p,2),clipped_panels(p,4)]
   af=plot(Xj,Yj,'LineWidth', 2.5)
   af.Color=[0 0 0];
end


end


