clear variables; close all;
%% INFO ===================================================================
% All orbits are considered in the same plane => (inclination and RAAN are the
% same for both orbits)
% 0. Define centeral planet/star/body of orbit
% 1. Define your target orbit (generally circular) eccentricity, semi-major
% axis and argument of perigee
% 2. Define your chaser orbit eccentricity, semi-major axis and argument of
% perigee
% 3. Choose the initial position of the target by defining the true anomaly
% 4. Choose the initial position of the chaser by defining the true anomaly
% 5. RUN the script
%==========================================================================
%% Centeral body ==========================================================
Cen_bod="SUN"; % or EARTH;

if Cen_bod=="EARTH"
    GM = 3.986005*10^5; % Earth gravitational parameter (km3/s2)
    R_cb=6378; % Earth radius [km]
elseif Cen_bod=="SUN"
    GM = 1.327*10^11; % Sun gravitational parameter (km3/s2)
    R_cb=695510; % Sun radius [km]
end

%% Target orbit parameters
e_t=0; % eccentricity
a_t=1.0816*10^8; % semi-major axis [km]
w_t=deg2rad(0); % argument of periapsis [°]

%% Chaser orbit parameters
e_c=0.02285; % eccentricity --> abs(1-(1-e_t)*a_t/(6378+300)) for dV burn at perigee
a_c=1.057433*10^8; % semi-major axis [km]
w_c=deg2rad(180); % argument of periapsis [°]

%% Initial Target configuration
theta_ti=deg2rad(60); %initial true anomaly [°]

%% Initial Chaser coniguration
theta_ci=deg2rad(-180); %initial true anomaly [°]

%% Calculate Target/Chaser position over time and LVLH

% Orbits of target and chaser
[X_t,Y_t]=plot_orbit(e_t,a_t,w_t,0,0,GM,1);
[X_c,Y_c]=plot_orbit(e_c,a_c,w_c,0,0,GM,1);

disp('Orbits of target and chaser calculated');

V=[];
H=[];
X=[];
Y=[];

T=0:100000:100000000; % time of plot in seconds

PHI_c=theta_ci+w_c;
PHI_t=theta_ti+w_t;

for t=T
    [x_t,y_t,phi_t]=plot_orbit(e_t,a_t,w_t,theta_ti,t,GM,2);
    [x_c,y_c,phi_c]=plot_orbit(e_c,a_c,w_c,theta_ci,t,GM,2);
    
    X=[X, [x_t;x_c]];
    Y=[Y, [y_t;y_c]];
    
    while PHI_t(end)-phi_t > 0.0001
        phi_t=phi_t+2*pi;
    end
    
    while PHI_c(end)-phi_c > 0.0001
        phi_c=phi_c+2*pi;
    end
    
    
    PHI_t=[PHI_t, phi_t];
    PHI_c=[PHI_c, phi_c];
    
    z_t=sqrt(x_t^2+y_t^2);
    z_c=sqrt(x_c^2+y_c^2);
    
    H = [H, arclength(phi_t,phi_c,a_t,a_t*sqrt(1-e_t^2))];
    
    V=[V,z_c-z_t];
    
end

%% Plot Target/Chaser position and LVLH over time

disp('Chaser is in blue on red orbit, while target is in red on blue orbit');

fig=figure('Position',[200 150 1000 500]);
set(fig,'defaultLegendAutoUpdate','off');
win(1) = subplot(1, 2, 1);
win(2) = subplot(1, 2, 2);

set(win,'Nextplot','add')

%Orbits
p1=plot(win(1),X_t,Y_t,'b');
p2=plot(win(1),X_c,Y_c,'r');


[xcb,ycb]=plotcircle(0,0,R_cb,0,360);
fill(win(1),xcb,ycb,'green') % Central body to scale

xL = xlim(win(1));
yL = ylim(win(1));
line(win(1),[0 0], yL,'Color','black','Linewidth',0.5);  %x-axis
line(win(1),xL , [0 0],'Color','black','Linewidth',0.5);  %y-axis

text(win(1),0,0,Cen_bod,'HorizontalAlignment','center','Color','black','FontSize',14)

axis(win(1),'image')
title(win(1),'Target and chaser orbits')
ylabel(win(1),'Y-ECI (km)');
xlabel(win(1),'X-ECI (km)');


%LVLH
plot(win(2),H,round(V),'r')
xL = xlim(win(2));
if 0>xL(2)
    xL(2)=0;
end
if 0<xL(1)
    xL(1)=0;
end

yL = ylim(win(2));
if 0>yL(2)
    yL(2)=0;
end
if 0<yL(1)
    yL(1)=0;
end

line(win(2),[0 0], yL,'Color','black','Linewidth',0.5);  %x-axis
line(win(2),xL , [0 0],'Color','black','Linewidth',0.5);  %y-axis

%axis(win(2),'image')
title(win(2),'LVLH (target centered)');
ylabel(win(2),'Vertical distance vs. target (km)');
xlabel(win(2),'Horizontal distance vs. target (km)');

%initialize target/chaser indicators
a1=plot(win(1),0,0);
a2=plot(win(1),0,0);
a3=plot(win(2),0,0);
t1=text(0,0,'');
t2=text(0,0,'');

i=0;
for t=T
    i=i+1;
    
    set(a1,'Visible','off');
    set(a2,'Visible','off');
    set(t2,'Visible','off');
    
    a1=plot(win(1),X(1,i),Y(1,i),'.r','MarkerSize',10,'DisplayName','Target');
    a2=plot(win(1),X(2,i),Y(2,i),'.b','MarkerSize',10);
    
    t2=text(win(1),0.01,0.02,string(t/(60*60*24))+' days','Units','normalized');
    
    set(a3,'Visible','off');
    set(t1,'Visible','off');
    
    a3=plot(win(2),H(i),round(V(i)),'.b','MarkerSize',10);
    
    t1=text(win(2),0.01,0.02,string(t/(60*60*24))+' days','Units','normalized');
    
    pause(0.01)
end
