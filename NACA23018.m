% function [cl,cm,control_x,control_y,cp,length]=NACA23018(n,a,v)%n为网格划分数，a为攻角（弧度制），v为无穷来流速度
close all;clc;clear all;
%     m=0.2025;k1=15.957;t=0.18;c=1;
n=100;a=0;v=10;
m=0.06;p=0.4;t=0.21;c=1;n_slot=10;n_trans=6;
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
trans_inj=[linspace(inj_point(1),x_u(sst_start),n_trans);linspace(inj_point(2),y_u(sst_start)-SST,n_trans)];
trans_suc=[linspace(x_u(sst_end),suc_point(1),n_trans);linspace(y_u(sst_end)-SST,suc_point(2),n_trans)];
%将上下翼面离散的点按顺时针方向合并
discrete_x=[x_u(1:inj_index),inj_points(1,2:n_slot),trans_inj(1,2:n_trans-1),x_u(sst_start:sst_end),trans_suc(1,2:n_trans-1),suc_points(1,:),x_u(suc_index+1:n),x_l(n-1:-1:1)];
discrete_y=[y_u(1:inj_index),inj_points(2,2:n_slot),trans_inj(2,2:n_trans-1),y_u(sst_start:sst_end)-SST,trans_suc(2,2:n_trans-1),suc_points(2,:),y_u(suc_index+1:n),y_l(n-1:-1:1)];
%     discrete_x=[x_u,x_l(n-1:-1:1)];
%     discrete_y=[y_u,y_l(n-1:-1:1)];
Flag=[ones(1,inj_index-1),2*ones(1,n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start),-2*ones(1,n_slot-1),ones(1,n-1+n-suc_index)];

%求控制面的切向与X轴夹角
mask_lamda=[ones(1,inj_index),zeros(1,2*n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start+1),zeros(1,2*n_slot-1),ones(1,n+n-suc_index)];
mask_face_type=[ones(1,inj_index-1),zeros(1,n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start),zeros(1,n_slot-1),ones(1,n-1+n-suc_index)];
mask_jet_type=[zeros(1,inj_index-1),ones(1,n_slot-1),zeros(1,2*n_trans-2+sst_end-sst_start),ones(1,n_slot-1),zeros(1,n-1+n-suc_index)];
figure(1);
plot(discrete_x,discrete_y);
hold on;
 
%     scatter(inj_point(1),inj_point(2),'*');
%     scatter(suc_point(1),suc_point(2),'*');
hold off;
xlim([0 1]);
ylim([-0.5 0.5]);
hold on;
control_x=zeros(1,length(discrete_x)-1);
control_y=zeros(1,length(discrete_y)-1);
Angle=zeros(1,length(Flag));
%控制点位置,从上翼面开始顺时针排列
for i=1:length(discrete_y)-1
    control_x(i)=(discrete_x(i+1)+discrete_x(i))/2;
    
    control_y(i)=(discrete_y(i+1)+discrete_y(i))/2;
 
    Angle(i)=atan2(discrete_y(i+1)-discrete_y(i),discrete_x(i+1)-discrete_x(i));
 
end

%对中弧线进行分段离散,计算每一段的长度
len=zeros(1,length(discrete_x)-1);
for i=1:length(discrete_x)-1
   len(i)=sqrt((discrete_x(i+1)-discrete_x(i))^2+ (discrete_y(i+1)-discrete_y(i))^2);
%    len(2*n-1-i)= sqrt((x_l(i+1)-x_l(i))^2+ (y_l(i+1)-y_l(i))^2);
end
V=[v*cos(a),v*sin(a)];
[results,B,D,E]=solver(V,control_x,control_y,discrete_x,discrete_y,n,Angle,Flag);



% 计算切向速度
% wing_suf_index=find(mask_lamda==1);
   
% wing_suf_vortex=lamda(wing_suf_index);
% wing_suf_vortex=(wing_suf_vortex(1:length(wing_suf_vortex)-1)+wing_suf_vortex(2:length(wing_suf_vortex)))/2;
panel_value=D*results;
vortex=panel_value.*(mask_face_type');
wing_suf_ind=find(mask_face_type==1);
%计算压力系数，并作图
figure(3);

 

vel=B*results+((V*[cos(Angle);sin(Angle)])').*mask_face_type'+((V*[cos(Angle+pi/2);sin(Angle+pi/2)])').*mask_jet_type';
wing_suf_x=discrete_x(wing_suf_ind);

check=E*results+((V*[cos(Angle);sin(Angle)])').*mask_jet_type'+((V*[cos(Angle+pi/2);sin(Angle+pi/2)])').*mask_face_type';
cp=(1-(vel.^2)/v^2);
plot(control_x(1:length(control_y)-n+1),cp(1:length(control_y)-n+1),'r.');
hold;
plot(control_x(length(control_y)-n:length(control_y)),cp(length(control_y)-n:length(control_y)),'b.');
set(gca,'YDir','reverse');
legend('上表面','下表面');
%计算升力系数
xlabel('弦长');
ylabel('压力系数');
title('压力系数曲线');
    

Gamma=len*vortex;%总环量
cl=2*Gamma/(v*c);
%计算力矩系数
% cm=-2*control_x.*len.*cos(Angle)*(lamda(1:2*n-2)+lamda(2:2*n-1))/2/(v*c);
% cm=2*control_x.*length.*cos(Angle)*cp/(v*c);
% flow(v,a,lamda,discrete_x,discrete_y,Angle,Flag)


% N_points=length(discrete_x)-1;
% V_u=zeros(n-1,1);
% for i=1:(n-1)
%     x_loc=(x_u(i)+x_u(i+1))/2;y_loc=(y_u(i)+y_u(i+1))/2;
%     
%     V_x=0;V_y=0;
%     index_col=1;
%     angle=atan2(y_u(i+1)-y_u(i),x_u(i+1)-x_u(i));
%     for j=1:N_points  
%         if Flag(j)==1
%             
%             [n_a,n_b,t_a,t_b]=linear_cof(x_loc,y_loc,discrete_x(j),discrete_y(j),discrete_x(j+1),discrete_y(j+1),Angle(j),angle); 
%             lamda=results(index_col:index_col+1)';
% %             V_x=V_x+lamda*[t_a;t_b];
% %             V_y=V_x+lamda*[n_a;n_b];    
%             V_u(i)=V_u(i)+lamda*[t_a;t_b];
%        
%         else 
%             [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(x_loc,y_loc,discrete_x(j),discrete_y(j),discrete_x(j+1),discrete_y(j+1),Angle(j),angle); 
%             gama=results(index_col:index_col+2)';
% %             V_x=V_x+gama*[t_a;t_b;t_c];
% %             V_y=V_x+gama*[n_a;n_b;n_c];       
%             V_u(i)=V_u(i)+gama*[t_a;t_b;t_c];
%         end
%         if Flag(j)~=Flag(mod(j,N_points)+1)
%             index_col=index_col+1;
%         end            
%         index_col=index_col+abs(Flag(j));
%     end
%     V_u(i)=V_u(i)+dot(V,[cos(angle),sin(angle)]);
% 
% end
    