%  Panel Code in MATLAB
%
%  Coded by L. sankar, April 1997
%
clear
close all
%% Create the panels and find the influsence co-efficients
U_inf  = 1  ;
% n_pan  = 64 ; % Number of panels to use
% panels = n_panel_circle(n_pan) ;  % Define the number of approximation panels
% I      =(zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ; % Initialise influence

% Create airfoil panels using jowkowski
aoa_degrees = 00 ;                                   % Angle of attack in degrees
panels      = jowkowski_function_5_0(aoa_degrees) ; % Create an airfoil in panels
n_pan       = length(panels);                       % Number of panels
% n_pan       = 5;
% panels      = n_panel_circle(n_pan) ;  % Define the number of approximation panels
I           = (zeros(n_pan,n_pan)) ; Phi_i=zeros(n_pan,1) ;   % Initialise influence


% Read Angle of Attack
%
alpha = aoa_degrees/pi*180;
%
% read number of points on the upper side of airfoil
%
nu = n_pan/2
%
% read number of points on the lower side of airfoil
%
nl = n_pan/2
%
% read Flag that states if this airfoil is symmetric
% if isym > 0 then airfoil is assumed symmetric
%
isym = 1
%
% Read a scaling factor
% The airfoil y- ordinates will be multiplied by this factor
%
factor= 3;

if(isym>0)
  nl = nu;
end
%
% Allocate storage for x and y
%
x = zeros(1,100);
y = zeros(1,100);
%
% Read the points on the upper surface
%
numel(nu)
   x = panels(:,1);
   y = panels(:,2);
   
%    x_append=x(1);
%    y_append=y(1);
%    
%    x=[x;x_append];
%    y=[y;y_append];
   
%
% If the airfoil is not symmetric, read lower side ordinates too..
%

%
% Plot the airfoil on window #1
%
plot(x,y);
n=nu+nl-2;
A=zeros(n+1,n+1);
ds=zeros(1,n);
pi=4. * atan(1.0);
%
% Assemble the Influence Coefficient Matrix A
%
 for i = 1:n
   t1= x(i+1)-x(i);
   t2 = y(i+1)-y(i);
   ds(i) = sqrt(t1*t1+t2*t2);
 end
 
for j = 1:n
 a(j,n+1) = 1.0;
 for i = 1:n
   if i == j
	 
     fprintf('This is Ds : %f \n',ds(i))
     fprintf('This is the log bit of Ds : %f \n',(log(0.5*ds(i)) - 1.0) )
     
     a(i,i) = ds(i)/(2.*pi) *(log(0.5*ds(i)) - 1.0);
   else
     xm1 = 0.5 * (x(j)+x(j+1));
     ym1 = 0.5 * (y(j)+y(j+1));
     dx  = (x(i+1)-x(i))/ds(i);
     dy  = (y(i+1)-y(i))/ds(i);
     t1  = x(i) - xm1;
     t2  = y(i) - ym1;
     t3  = x(i+1) - xm1;
     t7  = y(i+1) - ym1;
     t4  = t1 * dx + t2 * dy;
     t5  = t3 * dx + t7 * dy;
     t6  = t2 * dx - t1 * dy;
     t1  = t5 * log(t5*t5+t6*t6) - t4 * log(t4*t4+t6*t6);
     t2  = atan2(t6,t4)-atan2(t6,t5);
     a(j,i) = (0.5 * t1-t5+t4+t6*t2)/(2.*pi);
   end
 end
a(n+1,1) = 1.0;
a(n+1,n) = 1.0;
end
%
% Assemble the Right hand Side of the Matrix system
%
rhs=zeros(n+1,1);
alpha = alpha * pi /180;
xmid=zeros(n,1);
for i = 1:n
  xmid(i,1) = 0.5 * (x(i) + x(i+1));
  ymid = 0.5 * (y(i) + y(i+1));
  rhs(i,1) = ymid * cos(alpha) - xmid(i) * sin(alpha);
end
gamma = zeros(n+1,1);
%
% Solve the syetm of equations
% In MATLAB this is easy!
%
gamma = a\rhs;

I_matrix=a
V_inf_matrix=rhs;

cp=zeros(n,1);
cp1=zeros(n,1);
%
% Open a file to write x vs. Cp and the Loads
%
% Change the file name below, to open a new file every time 
%
fid=fopen('cp4.dat','w');
fprintf(fid,'   X          CP\n\n');
for i = 1:n
cp(i,1) = 1. - gamma(i) * gamma(i);
cp1(i,1) = - cp(i,1);
xa    = xmid(i,1);
cpa = cp(i,1);
%
% Write x and Cp to the file
%
% The xa- coordinate is the center points of panel 'i'
% Cpa is the Cp value at that point
%
fprintf(fid,'%10.4f %10.4f\n',xa,cpa);
end
%
% Open a new figure and plot x vs. Cp
%
figure(2);
plot(xmid,cp1);
%
% Compute Lift and Drag Coefficients
%
cy = 0.0;
cx = 0.0;
cm = 0.0;
% We assume that the airfoil has unit chord
% we assume that the leading edge is at i = nl;
for i=1:n
dx = x(i+1) - x(i);
dy = y(i+1) - y(i);
% xarm is the moment arem , equals distance from
% the center of the panel to quarter-chord.
xarm = 0.5 * (x(i+1)+x(i))-x(nl)-0.25;
cy = cy - cp(i,1) * dx;
cx = cx + cp(i,1) * dy;
cm = cm - cp(i,1) * dx * xarm;
end
%
% Print Lift and Drag coefficients on the screen
%
cl = cy * cos(alpha) - cx * sin(alpha)
cd = cy * sin(alpha) + cx * cos(alpha)
cm
%
% Write lift and Drag coefficients to a file
%
fprintf(fid,'  CL        CD  CM\n');
fprintf(fid,'%10.4f %10.4f %10.4f\n', cl,cd,cm);
fclose(fid);