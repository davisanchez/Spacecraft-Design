clear all
close all

%% Eclipses 
R_J=7.1398*10^4;
r_perijove=0.8*R_J;
r_apojove=5.5*R_J;
u_J=1.27*10^8;

lambda=asin(R_J/(2*r_perijove)) %[rad]
lambda=lambda*(360/(2*pi))
a_orbit=(r_perijove+r_apojove)/2;
T_orbit=sqrt(a_orbit^3/u_J);
T_eclipse=(2*lambda/360)*T_orbit; %[s]
T_eclipse=T_eclipse/3600


%% Mass budget
