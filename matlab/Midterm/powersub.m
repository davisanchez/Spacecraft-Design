clear all
close all

%David Sanchez del Rio

%% Constants
P_req=300; %[W] required power
R_J=7.1398*10^4;    %[km]
r_perijove=0.8*R_J+R_J; %[km]
r_apojove=5.5*R_J+R_J;  %[km]
u_J=1.27*10^8;
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



%% Psa calculation

P_sa = (k/(T_orbit-T_eclipse))*((P_d*(T_orbit-T_eclipse))/eta_d+(P_e*T_eclipse)/eta_e); 
        %[W]
        
        
%% Array dimensioning

solar_constant_e = 1361;    %[W/m2]
solar_constant_j = 50.4;%[W/m2]
%let's try multijunction cells
eff=0.274;
theta=25; %[deg] inclination angle------> to change
sat_life=8; %[years]
degradation_year=0.025/8;   %[percent]
Id=0.77;
P_o=solar_constant_j*eff;   %[W/m2] ideal solar output

%% BOL and EOL calculations
%demander si tout ce qui suit est juste
P_BOL=P_o*(Id*cos(theta*(2*pi/360)));   %[W/m2] input needed power
P_EOL=P_BOL*((1-degradation_year)^sat_life); %[W/m2] supposition for now 
A_sa=P_sa/P_EOL;    %[m2] area of solar pannels


%% Battery capacity
N=2;    %number of batteries
DoD=0.80;
n=0.9   %transmission efficiency

C_r=(T_eclipse*P_e)/(N*n*DoD)

%% Weight and power budget

%Solar arrays
sol_array_weight = 0.04*P_req %[kg]

%Batteries
battery_weight = C_r/45 %[kg], for NiH2
battery_weight=4*battery_weight %2 batteries, demander

%Power control Unit
pcu_weight = 0.02*P_req   %[kg]

%Regulators/converters
rc_weight = 0.025*P_req    %[kg]
rc_power = 0.2*P_req %[W] converted power
        
%Wiring
wiring_weight = 0.04*m_dry    %[kg]
wiring_power = 0.04*P_req  %[W]

mass_tot=sol_array_weight+battery_weight+pcu_weight+rc_weight%+wiring_weight

%% Area offset
m_wet=3625; %[kg]
s = 0.25*(m_wet)^(1/3); %[m] linear dimension 

L_a=1.5*s+0.5*sqrt(A_sa/2);  %[m] solar array area offset

%solar array moment of inertia [kg*m2]
%attention: we consider 2 squared solar panels
I_ax=(L_a^2+A_sa/12)*sol_array_weight; %perpendicular to array face
I_ay=(L_a^2+A_sa/24)*sol_array_weight; %%perpendicular to array axis
I_aa=(A_sa/24)*sol_array_weight; %about array axis