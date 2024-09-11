function [u_a,w_a,u_b,w_b,u_c,w_c] = Quadratic_source_cof(X,Y,x1,y1,x2,y2,angle,theta)
% 计算影响系数 
% 将控制面细分为n段求数值积分
%控制点j坐标
Xj=(x1+x2)/2;
Yj=(y1+y2)/2;
%控制点i在控制面j参考下下坐标，代入求得的影响系数公式
X_r=(X-Xj)*cos(angle)+(Y-Yj)*sin(angle);
Y_r=-(X-Xj)*sin(angle)+(Y-Yj)*cos(angle);
s=sqrt((x2-x1)^2+(y2-y1)^2);
r1=sqrt((X_r+s/2)^2+Y_r^2);
r2=sqrt((X_r-s/2)^2+Y_r^2);
% angle1=atan(Y_r/(X_r+s/2));
% angle2=atan(Y_r/(X_r-s/2));
angle1=atan2(Y_r,(X_r+s/2));
angle2=atan2(Y_r,(X_r-s/2));
theta_p=theta-angle;%控制点处的切向方向在面源坐标系下的角度
if(Y_r==0 && X_r==0)
    angle2=pi;
end
%计算三次面源在控制点处的速度，面源坐标系下，使用源单元
cof_xx_u=1/(2*pi)*(-X_r*s-2*log(r2/r1)*(X_r^2/2-Y_r^2/2)+2*X_r*Y_r*(angle2-angle1));
cof_xx_w=1/(2*pi)*(-(X_r^2-Y_r^2)*(angle1-angle2)+2*X_r*Y_r*log(r2/r1)-Y_r*(-s));
cof_x_u=1/(2*pi)*(X_r*log(r1/r2)-s+Y_r*(angle2-angle1));
cof_x_w=1/(2*pi)*(Y_r*log(r1/r2)+X_r*(angle2-angle1));
cof_u=1/(2*pi)*log(r1/r2);
cof_w=1/(2*pi)*(angle2-angle1);



up_a=cof_xx_u/(s^2/2)+cof_x_u/(-s);
wp_a=cof_xx_w/(s^2/2)+cof_x_w/(-s);
up_b=cof_xx_u/(-s^2/4)+cof_u;
wp_b=cof_xx_w/(-s^2/4)+cof_w;
up_c=cof_xx_u/(s^2/2)+cof_x_u/(s);
wp_c=cof_xx_w/(s^2/2)+cof_x_w/(s);
% cof=dot([up,wp],[cos(theta_p),sin(th
% eta_p)]);
u_a=up_a*cos(angle)-wp_a*sin(angle);
w_a=up_a*sin(angle)+wp_a*cos(angle);
u_b=up_b*cos(angle)-wp_b*sin(angle);
w_b=up_b*sin(angle)+wp_b*cos(angle);
u_c=up_c*cos(angle)-wp_c*sin(angle);
w_c=up_c*sin(angle)+wp_c*cos(angle);
%将面源坐标系下的速度投影到控制点的法向
end
