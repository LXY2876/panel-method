% function [cl,cm,control_x,control_y,cp,length]=NACA23018(n,a,v)%n为网格划分数，a为攻角（弧度制），v为无穷来流速度
    close all;clc;
    n=100;a=0;v=10;
%     m=0.2025;k1=15.957;t=0.18;c=1;
    m=0.06;p=0.4;t=0.21;c=1;
    syms X;
%     x1=linspace(0,m,n/2+1);
%     x2=linspace(m,c,n/2);
%     x=[x1(1:end-1),x2];
    x=linspace(0,c,n);
    
    
    %分段中弧线方程
%     Yc_f=k1/6*(X.^3-3*m*X.^2+m^2*(3-m).*X);
%     Yc_b=k1*m^3/6*(1-X);
    Yc_f=m*X/p^2*(2*p-X/c);
    Yc_b=m*(c-X)/(1-p)^2*(1-2*p+X/c);
    %机翼厚度
 
    y_t=t/0.2*c*(0.2969*(x./c).^0.5-0.126*x./c-0.3516*(x./c).^2+0.2843*(x./c).^3-0.1036*(x./c).^4);
   
    %对中弧线方程关于X求导
    dYc_f_dx=diff(Yc_f);
    dYc_b_dx=diff(Yc_b);
    
%     %根据网格数在中弧线上离散取点  
%     y_c=subs(Yc_f,X,x).*(x>=0 & x<m)+subs(Yc_b,X,x).*(x>=m & x<=1);
%     %根据网格数在中弧线导数上离散取点  ,并求出切线与X轴夹角，为弧度制
%     theta=double(atan(subs(dYc_f_dx,X,x).*(x>=0 & x<m)+subs(dYc_b_dx,X,x).*(x>=m & x<=1)));
    %根据网格数在中弧线上离散取点  
    y_c=subs(Yc_f,X,x).*(x>=0 & x<p*c)+subs(Yc_b,X,x).*(x>=p*c & x<=c);
    %根据网格数在中弧线导数上离散取点  ,并求出切线与X轴夹角，为弧度制
    theta=double(atan(subs(dYc_f_dx,X,x).*(x>=0 & x<p*c)+subs(dYc_b_dx,X,x).*(x>=p*c & x<=c)));
    %求出对应点的上下翼面位置
    y_u=double(y_c+y_t.*cos(theta));
    y_l=double(y_c-y_t.*cos(theta));
    x_u=double(x-y_t.*sin(theta));
    x_l=double(x-y_t.*sin(theta));
    %将上下翼面离散的点按顺时针方向合并
    discrete_x=[x_u,x_l(n-1:-1:1)];
    discrete_y=[y_u,y_l(n-1:-1:1)];
    %求控制面的切向与X轴夹角
    Angle=zeros(1,2*n-2);
    figure(1);
    plot(discrete_x,discrete_y);
    hold on;
    plot( x,y_c);
    hold off;
    xlim([0 1]);
    ylim([-0.2 0.2]);
    control_x=zeros(1,2*n-2);
    control_y=zeros(1,2*n-2);
    hold on;
    %控制点位置,从上翼面开始顺时针排列
    for i=1:2*n-2
        control_x(i)=(discrete_x(i+1)+discrete_x(i))/2;
%         control_x(2*n-1-i)=(x_l(i+1)+x_l(i))/2;
        control_y(i)=(discrete_y(i+1)+discrete_y(i))/2;
%         control_y(2*n-1-i)=(y_l(i+1)+y_l(i))/2;
        Angle(i)=atan2(discrete_y(i+1)-discrete_y(i),discrete_x(i+1)-discrete_x(i));
%         Angle(2*n-1-i)=atan2(y_l(i+1)-y_l(i),x_l(i+1)-x_l(i));
%         Angle(2*n-1-i)=atan2(y_l(i)-y_l(i+1),x_l(i)-x_l(i+1));
    end
    
    %对中弧线进行分段离散,计算每一段的长度
    length=zeros(1,2*n-2);
    for i=1:2*n-2
       length(i)=sqrt((discrete_x(i+1)-discrete_x(i))^2+ (discrete_y(i+1)-discrete_y(i))^2);
%        length(2*n-1-i)= sqrt((x_l(i+1)-x_l(i))^2+ (y_l(i+1)-y_l(i))^2);
    end
    V=[v*cos(a),v*sin(a)];
    [lamda,B]=solver(V,control_x,control_y,discrete_x,discrete_y,n,Angle);
     
    % 计算切向速度
     
   vel=B*lamda+(V*[cos(Angle);sin(Angle)])';
%     vel=(lamda(1:2*n-2)+lamda(2:2*n-1))/2;
    %计算压力系数，并作图
   figure(3);
   
   
   cp=(1-(vel.^2)/v^2);
   plot(control_x(1:n-1),cp(1:n-1),'r.');
   hold;
   plot(control_x(n:2*n-2),cp(n:2*n-2),'b.');
   set(gca,'YDir','reverse');
   legend('上表面','下表面');
   %计算升力系数
   xlabel('弦长');
   ylabel('压力系数');
   title('压力系数曲线');
        
   
   Gamma=length*(lamda(1:2*n-2)+lamda(2:2*n-1))/2;%总环量
   cl=2*Gamma/(v*c);
   %计算力矩系数
%    cm=-2*control_x.*length.*cos(Angle)*lamda/(v*c);
% cm=2*control_x.*length.*cos(Angle)*cp/(v*c);
   flow(n,v,a,lamda,control_x,control_y,discrete_x,discrete_y,Angle);
 
    