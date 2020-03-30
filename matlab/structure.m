clc
clear all
close all
%David Sanchez del Rio

%% Constants
g=9.81; %[m/s^2]
m_dry=1553.5; %[kg]
m_wet = 3825;   %[kg]
rho_al = 2.8*10^3; %[kg/m3] density of aluminium

%% Spacecraft size
spacecraft_volume = 0.025*m_wet %[m^3]
s = 0.25*(m_wet)^(1/3); %[m] linear dimension 
A_b = s^2;  %[m^2] body area
I_spacecraft = 0.01*(m_wet)^(5/3)  %[kg*m^2] moment of inertia
sat_height = 0.39*(m_wet)^(1/3); %[m]
sat_width = sqrt(A_b);    %[m]

buck = (10/4)*g*m_wet;   %[kg*m/s^2] max acceleration at launch
buck = 1.25*buck;   %25% margin
tank_vol = 1.1; %[m3]
tank_radius = (0.75*tank_vol/pi)^(1/3)  %[m]
tank_height = 4*tank_radius     %[m] total height for 2 tanks 

% aluminium 7075
E=71*10^9; %[N/m^2] young modulus

L = 3.0;    %[m] length of beams that will support the acceleration
L2=2;   %[m] length of beams that won't support the acceleration
I = buck*(L^2)/(pi^2*E)     %moment of inertia of beams
b=0.05; %we suppose beams are 5cm large
h = ((12*I)/b)^(1/3)    %[m] height of the beams

mass_structure = 4*b*h*L*rho_al + 4*b*h*L2*rho_al      %[kg], 8 beams

mass_plaques = (2*b+L2)*(L+2*h)*1*4 +(2*b+L2)*L2*1*2     %[kg] density of 1kg/m^2

total_mass_struc = mass_structure+mass_plaques  % 15% of m_dry ---> ok
