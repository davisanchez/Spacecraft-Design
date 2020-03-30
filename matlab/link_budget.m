%link budget
% David Sanchez del Rio

clear all; close all;

% Constants
re = 6378; % Earth radius km
k = 1.38064852*10^(-23); % Boltzmann constant m^2kg/s^2K

%% UP LINK

% Transmitting SGT
f_ul=7160*10^6; % frequency Hz
A = 591*10^6; % Satellite altitude km (Earth-Jupiter)
teta_e = deg2rad(50); % Elevation angle
alpha = asin(re*sin(pi/2+teta_e)/(re+A));
R = sqrt((re+A)^2+re^2-2*(re+A)*re*cos(pi/2-teta_e-alpha)) % Range km 
lambda = (299792458/f_ul);  % wavelength m 
Lp=10*log10((lambda/(4*pi*R*10^3))^2); % Free space path loss dB (negative)
DT= 34; % Antenna dish diameter m
GT=7.3466+20*log10(DT/lambda); % [dBi]
PT=10*log10(2000); % Power transmission [dBW]
EIRP=PT+GT; % Peak EIRP [dBW]
L = 0 ; % Uplink losses dB
G_1 = 10*log10(4*pi/lambda^2); % Gain of a 1m^2 antenna dB/m^2

OFD = EIRP + Lp + G_1 - L % Operating flux density dBW/m^2
SFD = -90 % Saturation flux density dBW/m^2
IBO = SFD-OFD % Input backoff dbW/m^2


if IBO < 0
    display('Satellite is being driven too hard')
else
    display('SFD > OFD, all is good')
end

% Receiving SAT 

DR = 2.4 ; 
T = 10*log10(250); % Noise temperature °K   à checker
GR = 7.3466+20*log10(DR/lambda); % Receive antenna gain [dBi]
G_T = GR-T; % Satellite G/T dB/K
C_T = EIRP + Lp - L + G_T ; % Satellite C/T dBW/K
C_N0 = C_T - 10*log10(k); 

%% DOWN LINK

% Transmitting SAT
f_dl=8412.3*10^6; % frequency Hz
lambda_dl = (299792458/f_dl);  % wavelength m 
PT_dl = 10*log10(50); % Transmission power [dBW] consider 50 W of output power
EIRP_dl = PT_dl+G_T; % EIRP transmit dBW

Lp_dl=10*log10((lambda_dl/(4*pi*R*10^3))^2); % Free space path loss dB (negative)

L_dl = 0; % Downlink losses dB

% Receiving SGT 
DR_dl = 70; 
T_dl=0;
GR_dl = 7.3466+20*log10(DR_dl/lambda_dl) % Receive antenna gain [dBi]
G_T_dl = GR_dl-T_dl % Satellite G/T dB/K    voir si mettre miscelanious 
C_T_dl = EIRP_dl + Lp_dl - L_dl + G_T_dl % Satellite C/T dBW/K
C_N0_dl = C_T_dl - 10*log10(k)

if C_N0_dl < C_N0
    display('System is downlink limited')
end

%Verify remaining margin is > 0
S_N_req = 4.4; % S/N required (Eb/N0 required) dB     voir cours pour la modulation
Lm = 6; % Miscellaneous losses (Implementation losses) dB
L_om = 3; % Operating margin for losses (Required link margin) dB
R_d = 3*10^3; % Data rate bits/s (bps)

S_N = C_N0_dl-10*log10(R_d)-S_N_req-Lm-L_om % Remaining margin dB

if S_N > 0
    display('Link is closed')
end






