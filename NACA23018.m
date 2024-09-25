% function [cl,cm,control_x,control_y,cp,length]=NACA23018(n,a,v)%n为网格划分数，a为攻角（弧度制），v为无穷来流速度
    clear ;close all;clc;
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
%     discrete_x=[inj_points(1,1:n_slot),trans_inj(1,2:n_trans-1),x_u(sst_start:sst_end),trans_suc(1,2:n_trans-1),suc_points(1,:),x_u(suc_index-1:-1:inj_index)];
%     discrete_y=[inj_points(2,1:n_slot),trans_inj(2,2:n_trans-1),y_u(sst_start:sst_end)-SST,trans_suc(2,2:n_trans-1),suc_points(2,:),y_u(suc_index-1:-1:inj_index)];
%     Flag=[2*ones(1,n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start),-2*ones(1,n_slot-1),-1*ones(1,suc_index-inj_index)];
    %直管道]\
    
%     discrete_x=[zeros(1,n_slot),linspace(0.1,9.9,2*n_trans-2+sst_end-sst_start),10*ones(1,n_slot),linspace(9.9,0.1,suc_index-inj_index),0];
%     discrete_y=[linspace(1,0,n_slot),linspace(-0.01,-0.49,2*n_trans-2+sst_end-sst_start),linspace(-0.5,1.5,n_slot),linspace(1.49,1.01,suc_index-inj_index),1];
    Flag=[2*ones(1,n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start+1),-2*ones(1,n_slot-1),-1*ones(1,suc_index-inj_index+1)];
    %弯管道
    discrete_x=[linspace(1.3,1,n_slot)*cos(2/3*pi),cos(linspace(11/18*pi,7/18*pi,2*n_trans-2+sst_end-sst_start)),linspace(1,1.3,n_slot)*cos(1/3*pi),1.3*cos(linspace(7/18*pi,12/18*pi,suc_index-inj_index+1))];
    discrete_y=[linspace(1.3,1,n_slot)*sin(2/3*pi),sin(linspace(11/18*pi,7/18*pi,2*n_trans-2+sst_end-sst_start)),linspace(1,1.3,n_slot)*sin(1/3*pi),1.3*sin(linspace(7/18*pi,12/18*pi,suc_index-inj_index+1))];
    mask_face_type=[zeros(1,n_slot-1),ones(1,2*n_trans-2+sst_end-sst_start),zeros(1,n_slot-1),ones(1,suc_index-inj_index)];
    mask_jet_type=[ones(1,n_slot-1),zeros(1,2*n_trans-2+sst_end-sst_start),ones(1,n_slot-1),zeros(1,suc_index-inj_index)];
    N_panels=length(Flag);
    %求控制面的切向与X轴夹角
    discrete_x1=discrete_x(1:N_panels);
    discrete_x2=discrete_x(2:N_panels+1);
    discrete_y1=discrete_y(1:N_panels);
    discrete_y2=discrete_y(2:N_panels+1);
    
    figure(1);
    plot(discrete_x,discrete_y);
    hold on;
 
    hold off;
    xlim([0 1]);
    ylim([-0.2 0.2]);
    control_x=(discrete_x1+discrete_x2)/2;
    control_y=(discrete_y1+discrete_y2)/2;
    Angle=atan2(discrete_y2-discrete_y1,discrete_x2-discrete_x1);
    hold on;
    %控制点位置,从上翼面开始顺时针排列
 
    
 
    V=[v*cos(a),v*sin(a)];
    [results,B,D,E]=solver(V,control_x,control_y,discrete_x,discrete_y,Angle,Flag);
    check=E*results;
    % 计算切向速度
    panel_value=D*results;
%    vel=B*lamda+(V*[cos(Angle);sin(Angle)])';
%     %计算压力系数，并作图
%    figure(3);
%    cp=(1-(vel.^2)/v^2);
%    plot(control_x(1:n-1),cp(1:n-1),'r.');
%    hold;
%    plot(control_x(n:2*n-2),cp(n:2*n-2),'b.');
%    set(gca,'YDir','reverse');
%    legend('上表面','下表面');
%    %计算升力系数
%    xlabel('弦长');
%    ylabel('压力系数');
%    title('压力系数曲线');
%    Gamma=length*(lamda(1:2*n-2)+lamda(2:2*n-1))/2;%总环量
%    cl=2*Gamma/(v*c);
   %计算力矩系数
%    cm=-2*control_x.*length.*cos(Angle)*lamda/(v*c);
% cm=2*control_x.*length.*cos(Angle)*cp/(v*c);
   flow(n,v,a,results,control_x,control_y,discrete_x,discrete_y,Angle,Flag);
 
    