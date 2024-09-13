function t_phi = linear_vortex_potential(x1,y1,x2,y2,angle,lamda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
syms x y
Xj=(x1+x2)/2;
Yj=(y1+y2)/2;
X_r=(x-Xj)*cos(angle)+(y-Yj)*sin(angle);
Y_r=-(x-Xj)*sin(angle)+(y-Yj)*cos(angle);
s=sqrt((x2-x1)^2+(y2-y1)^2);
 
angle1=atan2(Y_r,(X_r+s/2));
angle2=atan2(Y_r,(X_r-s/2));
r1=sqrt((X_r+s/2)^2+Y_r^2);
r2=sqrt((X_r-s/2)^2+Y_r^2);

phi_x=-1/(2*pi)*(X_r*Y_r*log(r1/r2)+Y_r/2*(-s)+(X_r^2-s^2/4-Y_r^2)*(angle1-angle2)/2);
phi=-1/(2*pi)*((X_r+s/2)*angle1-(X_r-s/2)*angle2+Y_r*log(r1/r2));
t_phi=lamda(1)*(phi_x/(-s)+phi/2)+lamda(2)*(phi_x/(s)+phi/2);

end