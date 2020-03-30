function L =arclength(t1,t2,a,b)



n=100;
tt = linspace(min(t1,t2),max(t1,t2),n);
dt = (max(t1,t2)-min(t1,t2))/(n-1);
dx = a*dt*sin(tt); 
dy = b*dt*cos(tt); 
L = sum(sqrt(dx.^2+dy.^2));

if t1<t2
    L=-L;
end

end
