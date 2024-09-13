function [normal_cof_a,normal_cof_b,tangl_cof_a,tangl_cof_b] = linear_cof(X,Y,x1,y1,x2,y2,angle,theta)
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
delta_theta=theta-angle;%控制点处的切向方向在面源坐标系下的角度
% cof=(X_r*log(r2/r1)+s*(log(r1*r2)-2)/2+Y_r*(angle2-angle1))/(2*pi);
%计算面源在控制点处的速度，面源坐标系下
%均匀面源
up_a=1/(2*pi)*(angle2-angle1)/2+1/(2*pi)*(Y_r*log(r2/r1)+X_r*(angle2-angle1))/(-s);
wp_a=1/(2*pi)*log(r2/r1)/2+1/(2*pi)*(X_r*log(r2/r1)+s-Y_r*(angle2-angle1))/(-s);
up_b=1/(2*pi)*(angle2-angle1)/2+1/(2*pi)*(Y_r*log(r2/r1)+X_r*(angle2-angle1))/(s);
wp_b=1/(2*pi)*log(r2/r1)/2+1/(2*pi)*(X_r*log(r2/r1)+s-Y_r*(angle2-angle1))/(s);
% cof=dot([up,wp],[cos(theta_p),sin(th
% eta_p)]);
if(Y_r==0 && X_r==0)
    up_a=0.25;
    wp_a=-1/(2*pi);
    up_b=0.25;
    wp_b=1/(2*pi);
end
% u_a=up_a*cos(angle)-wp_a*sin(angle);
% w_a=up_a*sin(angle)+wp_a*cos(angle);
% u_b=up_b*cos(angle)-wp_b*sin(angle);
% w_b=up_b*sin(angle)+wp_b*cos(angle);

%     u_a=up_a*cos(theta_p)-wp_a*sin(theta_p);
%     w_a=up_a*sin(theta_p)+wp_a*cos(theta_p);
normal_cof_a=dot([up_a,wp_a],[cos(delta_theta+pi/2),sin(delta_theta+pi/2)]);
%     u_b=up_b*cos(theta_p)-wp_b*sin(theta_p);
%     w_b=up_b*sin(theta_p)+wp_b*cos(theta_p);
normal_cof_b=dot([up_b,wp_b],[cos(delta_theta+pi/2),sin(delta_theta+pi/2)]);
tangl_cof_a=dot([up_a,wp_a],[cos(delta_theta),sin(delta_theta)]);
tangl_cof_b=dot([up_b,wp_b],[cos(delta_theta),sin(delta_theta)]);

%将面源坐标系下的速度投影到控制点的法向
end
