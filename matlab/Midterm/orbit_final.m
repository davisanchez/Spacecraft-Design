%David Sanchez del Rio
clear all
close all

%% constants and parameters ===============================================

%Earth
ue = 3.986*10^5; % Gravitational parameter [km^3/s^2]
Re = 6.378*10^3; % Equatorial radius [km]
ae = 1.496*10^8; % Semi-major axis of Earth orbit (heliocentric) [km]
g = 9.81;
%Sun
us = 1.32712442*10^11; % Gravitational parameter [km^3/s^2]
Ms = 1.9885*10^30; % Mass [kg]

aj=7.7832837*10^8; % Semi-major axis of Jupiter orbit (heliocentric) [km]

%% From the Earth to Jupiter

a=(ae+aj)/2; %semi major axis of orbital transfer [km]

V_d=sqrt(us*(2/ae-1/a)); %Heliocentric Earth departure velocity  [km/s]
V_e=sqrt(us/ae); % Earth heliocentric velocity [km/s]

V_d_inf=V_d-V_e; % Earth departure excess velocity [km/s]

P=2*pi*sqrt(a^3/us); % Elliptic orbital period of Hohmann transfer [s]
T=(P/2)/(3600*24*365); % Transfer time [years]



%% Jupiter position at Probe departure 
V_j=sqrt(us/aj); % Jupiter heliocentric velocity
alpha_J=180-360/((2*pi*aj)/V_j)*3600*24*365*T; % angle between Earth and Jupiter at launch [°]

%% Arrival velocity at infinity around Jupiter
V_a=sqrt(us*(2/aj-1/a)); % [km/s] Arrival heliocentric velocity
V_a_inf=V_j-V_a; %[km/s] Arrival excess velocity

%% Orbit, final delta_v

R_J=7.1398*10^4;    %[km]
r_perijove=0.8*R_J+R_J; %[km]
r_apojove=5.5*R_J+R_J;  %[km]
u_J=1.27*10^8;
a_J_ell=(r_perijove+r_apojove)/2;   %[km]

c3=V_d_inf^2;   %[km2/m2]
v_hyp_perijove=sqrt(V_a_inf^2+2*(u_J/r_perijove)); %[km/s]
v_ell_final=sqrt(u_J*(2/r_perijove-1/r_apojove)); %[km/s] at perijove
dv_orbit_insertion=v_hyp_perijove-v_ell_final;  % [km/s] breaking

delta_v_final=dv_orbit_insertion; %[km/s]

%% Change of plane (equatorial ot polar) at apojove
%let's suppose a 85 degrees plane change
launch_angle=28;
alpha=85-launch_angle;
e=0.5663;

v_apojove=sqrt(u_J*((1-e)/(r_apojove*1000)));  %[km/s]
dv_change_plane=2*v_apojove*sin((alpha*2*pi/360)/2); %[km/s]


%% Eclipses


lambda=asin(R_J/(r_perijove)); %[rad]
lambda=lambda*(360/(2*pi)) %[deg]
T_orbit=2*pi*sqrt(r_apojove^3/u_J); %[s]
T_eclipse=(2*lambda/360)*T_orbit; %[s]
T_eclipse=T_eclipse/3600;   %[h]


%% DeltaV budget

Isp_main=305;    %[s]
Isp_corr=240;  %[s]
m_wet=3625; %[kg]


corr = [139 , 50, 1081 , 27 , 747]; %[m/s], corrections
corr=1.1*corr %taking 10% margin

m_it = zeros (1 , length (corr) ) ;
%delta_T_it = zeros (1 , length (corr) ) ;

for i = 1: length (corr)
    if i == 3
        m_it( i ) = m_wet * (1 - exp(-corr( i ) / (g*Isp_main)) )  
        m_wet = m_wet - m_it( i ) ;
    else
        m_it( i ) = m_wet * (1 - exp(-corr( i ) / (g*Isp_corr)) ) ; 
        m_wet = m_wet - m_it( i ) ;
    end
end
m_wet=3625;
m_prop=sum(m_it);


% delta_v_final=2.044; %taking Juno's fredicted best delta v
% delta_v_final_main=2.044-0.2;   %main burn deltaV
% delta_v_final_corr1=0.05; %corrections deltaV
% delta_v_final_corr2=0.15; %corrections deltaV
% 
% delta_v_final=delta_v_final+0.1*delta_v_final;   %taking 10% margin
% delta_v_final_main=delta_v_final_main+0.1*delta_v_final_main;
% delta_v_final_corr1=delta_v_final_corr1+0.1*delta_v_final_corr1;
% delta_v_final_corr2=delta_v_final_corr2+0.1*delta_v_final_corr2;
% 
% %m_prop=m_wet*(-exp(-delta_v_final*1000/(g*Isp))+1);
% %m_prop=m_wet*(-exp(-(dv_change_plane+0.1*dv_change_plane)*1000/(g*Isp))+1);
% m_prop_corr1=m_wet*(-exp(-delta_v_final_corr1*1000/(g*Isp_corr))+1);
% m_prop_main=(m_wet-m_prop_corr1)*(-exp(-delta_v_final_main*1000/(g*Isp_main))+1);
% 
% m_prop_corr2=(m_wet-m_prop_main-m_prop_corr1)*(-exp(-delta_v_final_corr2*1000/(g*Isp_corr))+1);
% m_prop=m_prop_main+m_prop_corr1+m_prop_corr2; %total mass from main burn and correction burns
% 

mf=m_wet-m_prop;    %[kg] based on Juno's weight
F=mf*0.1*g; %thrust [N]
mass_propsyst=26.4 + 0.077*m_prop;  %[kg]
mass_propsyst=mass_propsyst+0.25*mass_propsyst; %taking worst case mass

%% Duration of burns
m_it = zeros (1 , length (corr) ) ;
delta_T_it = zeros (1 , length (corr) ) ;

for i = 1: length (corr)
    if i == 3 % the thirs item is the aerobraking (main) burn
        m_it( i ) = m_wet * (1 - exp(-corr( i ) / (g*Isp_main)) )  
        m_wet = m_wet - m_it( i ) ;
        delta_T_it(i)=((m_wet*g*Isp_main)/F)*(1-exp(-corr(i)/(g*Isp_main)));
    else
        m_it( i ) = m_wet * (1 - exp(-corr( i ) / (g*Isp_corr)) ) ; 
        m_wet = m_wet - m_it( i ) ;
        delta_T_it(i)=((m_wet*g*Isp_corr)/F)*(1-exp(-corr(i)/(g*Isp_corr)));
    end
end
% 
% delta_T_corr1=((m_wet*g*Isp_corr)/F)*(1-exp(-delta_v_final_corr1*1000/(g*Isp_corr))); %[s]
% delta_T_joi=(((m_wet-m_prop_corr1)*g*Isp_main)/F)*(1-exp(-delta_v_final_main*1000/(g*Isp_main))); %[s]
% delta_T_corr2=(((m_wet-m_prop_corr1-m_prop_main)*g*Isp_corr)/F)*(1-exp(-delta_v_final_corr2*1000/(g*Isp_corr))); %[s]
% %delta_T_=((m_wet*g*Isp)/F)*(1-exp(-delta_v_final*1000/(g*Isp))); %[s]
% delta_T=delta_T_joi+delta_T_corr1+delta_T_corr2; %[s]

%% Jupiter's coverage

epsilon=20;  %[deg] elevation angle
eta=30;  %[deg] instant coverage of the measurement devices--> to determine
lambda=90-(eta+epsilon); %[deg]
D=R_J*(sin(lambda*(2*pi/360))/sin(eta*(2*pi/360)));  %[km]
W_F=R_J*asin(D*sin(eta*(2*pi/360))/R_J);   %[km]

F_A=D*(pi/4)*W_F;   %instntaneous footprint area [km^2]
%F_A=(pi/4)*D^2*((sin(eta*2*pi/360))^2/sin(epsilon*2*pi/360))


%% Spacecraft size
spacecraft_volume = 0.01*m_wet; %[m^3]
s = 0.25*(m_wet)^(1/3); %[m] linear dimension 
A_b = s^2;  %[m^2] body area
I_spacecraft = 0.01*(m_wet)^(5/3);  %[kg*m^2] moment of inertia

%% Orbital characteristics of planets
% Earth
Trev_e=1; % Revolution period [years]
w_e = 360/Trev_e; % Angular speed in heliocentric rf [°/year]
teta_init_e = 10*30.5*(360/365); % Angle of Earth with respect to vernal equinox on 08.10.2019 at 00:00 UT [Â°]

% Jupiter
Trev_j=1/(V_j/aj)/(3600*24*365)*2*pi; % Revolution period [years]
w_j = 360/Trev_j; % Angular speed in heliocentric rf [°/year]
teta_init_j = 10*30.5*(360/(Trev_j*365)); % Angle of Jupiter with respect to vernal equinox on 01.10.2018 at 00:00 UT [Â°]

%% Calculating Launch date and time =======================================

Tsyn=abs((Trev_e*Trev_j)/(Trev_e-Trev_j)); % Synodic period [years]


alpha = 360*Trev_e/Trev_j; % Angle destination planet moves during travel time [°]
beta = abs(180-alpha); % Angle between departure and destination planet at launch [Â°]
beta_init = teta_init_e-teta_init_j; % Angle between earth and Jupiter on 01.10.2018 at 00:00 UT [Â°]

%T_launch=datetime('7-August-2011 00:00:00')+(beta-beta_init)/(w_e-w_j)*365;
T_launch=datetime('20-January-2025 00:00:00')+abs(beta-beta_init)/(w_e-w_j)*365;


for i=1:2
    T_launch=[T_launch, T_launch(end)+Tsyn*365]
end


