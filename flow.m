function [] = flow(v,a,results,discrete_x,discrete_y,theta,Flag)
%流场粒子仿真
    syms x y
    V=[v*cos(a),v*sin(a)];
%     [discrete_x,discrete_y]=angleChange(discrete_x,discrete_y,a);
%     [control_x,control_y]=angleChange(control_x,control_y,a);	
    % 流线细分设置
    div_y=0.07;
    % 流场模拟box边界
    x_l=-0.3;x_r=1.3;
    y_l=-0.5;y_u=0.5;
    % 时间步进分度
    dt=0.003;
    % y轴向偏移避开穿模
    y_shift=0.03;

    % 创建粒子发射源
    N=round((y_u-y_l)/div_y+1);
    y_l=y_l+y_shift;y_u=y_u+y_shift;
    x_emit=zeros(1,N)+x_l;
    y_emit=y_u:-1*div_y:y_l;

    % 解算流线坐标信息
    X=x_emit';Y=y_emit';
    i=1;
    N_points=length(discrete_x)-1;
    V_x=0;
    V_y=0;
    index_col=1;
    for j=1:N_points  
        if Flag(j)==1
            [n_a,n_b,t_a,t_b]=linear_cof(x,y,discrete_x(j),discrete_y(j),discrete_x(j+1),discrete_y(j+1),theta(j),0); 
            lamda=results(index_col:index_col+1)';
            V_x=V_x+lamda*[t_a;t_b];
            V_y=V_x+lamda*[n_a;n_b];
               
        else 
            gama=results(index_col:index_col+2)';
            [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(x,y,discrete_x(j),discrete_y(j),discrete_x(j+1),discrete_y(j+1),theta(j),0); 
            V_x=V_x+gama*[t_a;t_b;t_c];
            V_y=V_x+gama*[n_a;n_b;n_c];     
            
     
        end
        if Flag(j)~=Flag(mod(j,N_points)+1)
            index_col=index_col+1;
        end            
        index_col=index_col+abs(Flag(j));
    end

    V_x=simplify(V_x);
    V_y=simplify(V_y);
%     for k=1:N_points  
%     
%         if Flag(k)==1
%             lamda=results(index_col:index_col+1);
%             t_phi=t_phi+linear_vortex_potential(discrete_x(k),discrete_y(k),discrete_x(k+1),discrete_y(k+1),theta(k),lamda); 
%        
%         else 
%             gama=results(index_col:index_col+2);
%             t_phi=t_phi+Quadratic_source_potential(discrete_x(k),discrete_y(k),discrete_x(k+1),discrete_y(k+1),theta(k),gama); 
%   
%         end
%         if Flag(k)~=Flag(mod(k,N_points)+1)
%             index_col=index_col+1;
%         end            
%         index_col=index_col+abs(Flag(k));
%     end
%     V_x=diff(t_phi,x);
%     V_y=diff(t_phi,y);
    

    while min(X(:,i))<x_r && i<50
        for j=1:N
            x_loc=X(j,i);y_loc=Y(j,i);      
 


%             % 控制点到（x，y）的距离
%             dist2=abs((-control_x+x).^2+(-control_y+y).^2);
%             % 诱导速度矢量
%             v_x0=lamda'.*length.*(-control_y+y)./(2*pi*dist2);
%             v_y0=lamda'.*length.*(control_x-x)./(2*pi*dist2);
            % 面涡对（x，y）点产生的诱导速度v
            v_x=double(subs(V_x,{x,y},{x_loc,y_loc}))+V(1);
            v_y=double(subs(V_y,{x,y},{x_loc,y_loc}))+V(2);
            %下一时刻例子的位置
            X(j,i+1)=X(j,i)+v_x*dt; 
            Y(j,i+1)=Y(j,i)+v_y*dt;
        end 
        i=i+1;
    end

    % 粒子运动可视化


    index=1;
     
    [~,frame]=size(X);%frame为总帧数
    step=50;
    delta_t=0.01;%用于调整播放速度
    while index<frame
        hold off;
        figure(4);
        plot(discrete_x,discrete_y)
        axis equal
        axis([x_l x_r y_l y_u])
        hold on
        plot(X(:,mod(index,7)+1:7:end) ,Y(:,mod(index,7)+1:7:end) ,'b.','MarkerSize',13,'color',[77/255 190/255 238/255])
        plot(X(:,mod(index,step)+1:step:index+step),Y(:,mod(index,step)+1:step:index+step),'.','MarkerSize',13,'color',[217/255 83/255 25/255])
        index=index+1;
        pause(delta_t);
    end
end