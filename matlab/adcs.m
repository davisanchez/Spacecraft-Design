clear all;
close all;

%David Sanchez del Rio

%% Constants
P_req=300; %[W] required power
R_J=7.1398*10^4;    %[km]
r_perijove=0.8*R_J+R_J; %[km]
r_apojove=5.5*R_J+R_J;  %[km]
u_J=1.27*10^8;
ue = 3.986*10^5;
Re = 6.378*10^3;
k=1; %product of all degradations
eta_d=0.85; %solar array to loads efficiency day
eta_e=0.65; %solar array to loads efficiency eclipse
P_d = P_req; %[W] power day
P_e = P_req-30; %[W] power eclipse
m_dry=1553.5   %[kg]

lambda=asin(R_J/(r_perijove)); %[rad]
lambda=lambda*(360/(2*pi)) %[deg]
T_orbit=2*pi*sqrt(r_apojove^3/u_J); %[s]
T_eclipse=(2*lambda/360)*T_orbit; %[s]
T_eclipse=T_eclipse/3600;   %[h]
T_orbit=T_orbit/3600;   %[h]



%% Calculation of external torques
m_wet=3625; %[kg]
s = 0.25*(m_wet)^(1/3); %[m] linear dimension 
b=1;    %[m] from CoP to panels ----> à revoir
A_sr=47.75; %[m2]
P_sr = 3.4*10^(-7); %[N/m2]
C_r=0.9;
M_J=2.83*10^20; %[kg]
r=r_perijove;   %[m]
mom_sat = 1  

F_sr=A_sr*P_sr*(1+C_r)  %[N]
T_sr=b*F_sr     %[Nm]
B=2*M_J/((1000*r)^3)    %[T]
T_mag = mom_sat*B   %[Nm]

%atmospheric drag flyby
C_d=2.2;
A=3*2;
rho=4*10^-13;
vrel=sqrt(ue/(Re));

F_drag=0.5*C_d*A*rho*vrel^2;
tdrag=1.5*F_drag

%% Gravity gradient
y_height = 3;   %[m]
z_height = 2;   %[m]
x_height = 2;   %[m]
m_current = 2230.3; % mass after JOI
theta=1; %[deg]
theta = theta*(2*pi/360);   %[rad]
I_y_panel = 831.48; %[kg*m^2]
I_z_panel = 855.36; %[kg*m^2]
I_y_sat = (m_current/12)*(y_height^2+z_height^2); %[kg*m^2]
I_z_sat = (m_current/12)*(x_height^2+z_height^2); %[kg*m^2]
I_y = I_y_panel+I_y_sat
I_z=I_z_panel+I_z_sat

T_g=((3*1000*u_J)/(2*(1000*r_perijove)^3))*sin(2*theta)*abs(I_y_panel+I_y_sat-I_z_panel-I_z_sat)

dist = [T_g T_sr T_mag tdrag];

T_max = max(dist) %taking the maximal torque

%% Reaction wheels

h_reac = T_max*0.707*(T_orbit*3600)/4    %[Nms]
h_reac2 = T_max*0.707/4

%% Momentum wheels
theta_a=0.06/2;

h_mom = T_max*(T_orbit*3600)/(4*theta_a)    %[Nms]
h_mom2 = T_max/(4*theta_a)

%% Propellant mass
L=1.5; %[m] moment arm length 
Isp=240;
g=9.81;
t_pulse=1; %aucune idée
L_days = 3*365; %time post orbit in days
N_wheels=3;

F=(h_mom+h_reac)/L;

M_prop = (N_wheels*t_pulse*L_days*F)/(g*Isp)    %aucune idée encore




