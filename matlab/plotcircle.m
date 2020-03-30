function [xc,yc]=plotcircle(x,y,r,ang_i,ang_f)
%x and y are the coordinates of the center of the circle
%r is the radius of the circle
%0.01 is the angle step, bigger values will draw the circle faster but
%you might notice imperfections (not very smooth)
ang=ang_i:0.1:ang_f; 
xc=r*cosd(ang)+x;
yc=r*sind(ang)+y;
end