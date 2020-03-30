%thermal control
%David Sanchez del Rio

clc; clear all; close all;


Psun = 3.856*10^26; %[W] %solar rad pressure
d_SJ =  778300000; %[km]
d_SE = 149.6*10^6; %[km]
alpha=0.25;
epsilon=0.04;
surface=6;  %m^2
sigma = 5.67*10^(-8);
Arad = 4*6+2*4;

J_SJ = Psun/(4*pi*(1000*d_SJ)^2)

% Gravity assist
J_SE = Psun/(4*pi*(1000*d_SE)^2)
Q_solar_GA=alpha*surface*J_SE
J_a_GA=0.37; %albedo, SMAD
Q_al_GA=alpha*surface*J_a_GA
Q_IR_GA=231; %[W/m^2] 

%Orbit 
J_SJ = Psun/(4*pi*(1000*d_SJ)^2)
Q_solar_SJ=alpha*surface*J_SJ
J_a_SJ=0.34; %albedo, SMAD
Q_al_SJ=alpha*surface*J_a_SJ
Q_IR_SJ=13.5; %[W/m^2]

%black 
T_hot_black = ((Q_solar_GA+Q_al_GA+Q_IR_GA)/(epsilon*sigma*Arad))^(1/4)
T_cold_black = ((Q_IR_GA+Q_al_GA)/(epsilon*sigma*Arad))^(1/4)

%white 
T_hot_gold = ((Q_solar_GA+Q_al_GA+Q_IR_GA)/(epsilon*sigma*Arad))^(1/4)
T_cold_gold = ((Q_IR_GA+Q_al_GA)/(epsilon*sigma*Arad))^(1/4)

%gold 
T_hot_gold = ((Q_solar_GA+Q_al_GA+Q_IR_GA)/(epsilon*sigma*Arad))^(1/4)
T_cold_gold = ((Q_IR_GA+Q_al_GA)/(epsilon*sigma*Arad))^(1/4)


T_eq_GA = ((Q_solar_GA+Q_al_GA+Q_IR_GA)/(epsilon*sigma*Arad))^(1/4)
T_eq_SJ = ((Q_solar_SJ+Q_al_SJ+Q_IR_SJ)/(epsilon*sigma*Arad))^(1/4)