% function [cl,cm,control_x,control_y,cp,length]=NACA23018(n,a,v)%n为网格划分数，a为攻角（弧度制），v为无穷来流速度
close all;clc;
%     m=0.2025;k1=15.957;t=0.18;c=1;
n=100;a=0;v=10;
m=0.06;p=0.4;t=0.21;c=1;n_slot=10;
inj_size=0.01*c;suc_size=0.02*c;
inj_loc=0.07;suc_loc=0.78;SST=0.01*c;
inj_angle=40/180*pi;suc_angle=70/180*pi;
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

%根据网格数在中弧线上离散取点  
y_c=subs(Yc_f,X,x).*(x>=0 & x<p*c)+subs(Yc_b,X,x).*(x>=p*c & x<=c);
%根据网格数在中弧线导数上离散取点  ,并求出切线与X轴夹角，为弧度制
theta=double(atan(subs(dYc_f_dx,X,x).*(x>=0 & x<p*c)+subs(dYc_b_dx,X,x).*(x>=p*c & x<=c)));

%求出对应点的上下翼面位置
y_u=double(y_c+y_t.*cos(theta));
y_l=double(y_c-y_t.*cos(theta));
x_u=double(x-y_t.*sin(theta));
x_l=double(x-y_t.*sin(theta));
inj_index=round(inj_loc*n);
suc_index=round(suc_loc*n);

x_u_SST=x_u(inj_index:suc_index);
y_u_SST=y_u(inj_index:suc_index)-SST;
inj_point=[x_u(inj_index)+sin(inj_angle)*inj_size,y_u(inj_index)-cos(inj_angle)*inj_size];
suc_point=[x_u(suc_index)-sin(suc_angle)*suc_size,y_u(suc_index)-cos(inj_angle)*suc_size];
i=inj_index;
while inj_point(1)>=x_u(i)
    i=i+1;
end
sst_start=i+5;
i=suc_index;
while suc_point(1)<=x_u(i)
    i=i-1;
end
sst_end=i-5;
%吹气口和吸气口处加密
inj_points=[linspace(x_u(inj_index),inj_point(1),n_slot);linspace(y_u(inj_index),inj_point(2),n_slot)];
suc_points=[linspace(suc_point(1),x_u(suc_index),n_slot);linspace(suc_point(2),y_u(suc_index),n_slot)];
%过渡段
trans_inj=[linspace(inj_point(1),x_u(sst_start),n_slot);linspace(inj_point(2),y_u(sst_start)-SST,n_slot)];
trans_suc=[linspace(x_u(sst_end),suc_point(1),n_slot);linspace(y_u(sst_end)-SST,suc_point(2),n_slot)];
%将上下翼面离散的点按顺时针方向合并
discrete_x=[x_u(1:inj_index),inj_points(1,2:n_slot),trans_inj(1,2:n_slot-1),x_u(sst_start:sst_end),trans_suc(1,2:n_slot-1),suc_points(1,:),x_u(suc_index+1:n),x_l(n-1:-1:1)];
discrete_y=[y_u(1:inj_index),inj_points(2,2:n_slot),trans_inj(2,2:n_slot-1),y_u(sst_start:sst_end)-SST,trans_suc(2,2:n_slot-1),suc_points(2,:),y_u(suc_index+1:n),y_l(n-1:-1:1)];
%     discrete_x=[x_u,x_l(n-1:-1:1)];
%     discrete_y=[y_u,y_l(n-1:-1:1)];
Flag=[ones(1,inj_index-1),-1,2*ones(1,n_slot-1),-1,ones(1,2*n_slot-2+sst_end-sst_start),-1,-2*ones(1,n_slot-1),-1,ones(1,n-1+n-suc_index)];
%求控制面的切向与X轴夹角

figure(1);
plot(discrete_x,discrete_y);
hold on;
 
%     scatter(inj_point(1),inj_point(2),'*');
%     scatter(suc_point(1),suc_point(2),'*');
hold off;
xlim([0 1]);
ylim([-0.2 0.2]);
hold on;
control_x=zeros(1,length(discrete_x)-1);
control_y=zeros(1,length(discrete_y)-1);
Angle=zeros(1,length(Flag));
%控制点位置,从上翼面开始顺时针排列
for i=1:length(discrete_y)-1
    control_x(i)=(discrete_x(i+1)+discrete_x(i))/2;
    
    control_y(i)=(discrete_y(i+1)+discrete_y(i))/2;
 
    Angle(i)=atan2(discrete_y(i+1)-discrete_y(i),discrete_x(i+1)-discrete_x(i));
%         Angle(2*n-1-i)=atan2(y_l(i+1)-y_l(i),x_l(i+1)-x_l(i));
%     Angle(2*n-1-i)=atan2(y_l(i)-y_l(i+1),x_l(i)-x_l(i+1));
end

%对中弧线进行分段离散,计算每一段的长度
len=zeros(1,2*n-2);
for i=1:n-1
   len(i)=sqrt((x_u(i+1)-x_u(i))^2+ (y_u(i+1)-y_u(i))^2);
   len(2*n-1-i)= sqrt((x_l(i+1)-x_l(i))^2+ (y_l(i+1)-y_l(i))^2);
end
V=[v*cos(a),v*sin(a)];
[lamda,B]=solver(V,control_x,control_y,discrete_x,discrete_y,n,Angle,Flag);
 
% 计算切向速度
wing_index=find(Flag==1);
%    vel=B*lamda+(V*[cos(Angle);sin(Angle)])';
wing_suf_vel=(lamda(wing_index)+lamda(wing_index+1))/2;
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
    

Gamma=len*(lamda(1:2*n-2)+lamda(2:2*n-1))/2;%总环量
cl=2*Gamma/(v*c);
%计算力矩系数
cm=-2*control_x.*len.*cos(Angle)*(lamda(1:2*n-2)+lamda(2:2*n-1))/2/(v*c);
% cm=2*control_x.*length.*cos(Angle)*cp/(v*c);
%    flow(n,v,a,lamda,control_x,control_y,discrete_x,discrete_y,length);
   
   
    