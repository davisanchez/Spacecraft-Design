function [X,Y,phi]=plot_orbit(e,a,w,thetai,t_o,GM,arg)

T=2*pi*sqrt(a^3/GM);
n=2*pi/T; %mean motion

Ei=acos((e+cos(thetai))/(1+e*cos(thetai)));

ti=(Ei-e*sin(Ei))/n;
if mod(thetai,2*pi)>=pi
ti=T-ti;
end


X=[];
Y=[];

if arg==1
    Tend=T;
    m=1000;
else
    Tend=ti+t_o;
    m=1;
end

for t=linspace(ti+t_o,Tend,m)

% Mean anomaly at ti
M =(sqrt(GM/a^3))*t;

% Eccentric anomaly

E1=M; % initialize E
error=1; % initialize error
while error>=10^-3
    E2=M+e*sin(E1);
    error=abs(E2-E1);
    E1=E2;
end

E=E2;

% True anomaly
f = atan2((sqrt(1-e^2)*sin(E)),(cos(E)-e));

% Argument of perigee with true anomaly
phi = w+f;

% Radial distance
r = a*(1-e.*cos(E));

% Cartesian coordinates

x = r*cos(phi);
y = r*sin(phi);

X=[X,x]; % Make vector of x coordinates
Y=[Y,y]; % Make vector of y coordinates
end

end
