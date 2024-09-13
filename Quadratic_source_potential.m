function t_phi = Quadratic_source_potential(x1,y1,x2,y2,angle,sigma)
%UNTITLED2 Summary of this function goes here
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
phi_xx=1/(2*pi)*(2*Y_r^2*s/3-2*Y_r^3*log(r2/r1)/3+2*((s/2)^3*(log(r2)+log(r1)))-1/18*s^3+2*X_r*Y_r^2*log(r2/r1)+2*Y_r^3*(angle1-angle2)/3+2*X_r^2*Y_r*(angle2-angle1));
phi_x=1/(4*pi)*((X_r^2-s^2/2-Y_r^2)*log(r1/r2)+2*X_r*Y_r*(angle2-angle1)-X_r*s);
phi=1/(2*pi)*((X_r+s/2)*log(r1)-(X_r-s/2)*log(r2)+Y_r*(angle2-angle1));
t_phi=sigma(1)*(phi_xx/(s^2/2)+phi_x/(-s))+sigma(2)*(phi_xx/(-s^2/4)+phi)+sigma(3)*(phi_xx/(s^2/2)+phi_x/(-s));

end