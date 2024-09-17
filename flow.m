function [] = flow(n,v,a,results,control_x,control_y,discrete_x,discrete_y,theta)
%流场粒子仿真
    V=[v*cos(a),v*sin(a)];
%     [discrete_x,discrete_y]=angleChange(discrete_x,discrete_y,a);
    [control_x,control_y]=angleChange(control_x,control_y,a);	
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
    x_emit=zeros(1,N)+0.2;
    y_emit=[y_u:-1*div_y:y_l];

    % 解算流线坐标信息
    X=x_emit';Y=y_emit';
    k=1;
    while min(X(:,k))<x_r  && k<200
        for i=1:N
            x_loc=X(i,k);y_loc=Y(i,k);    
            V_x=V(1);V_y=V(2);
             for j=1:2*n-2
           
                [n_a,n_b,t_a,t_b]=linear_cof(x_loc,y_loc,discrete_x(j),discrete_y(j),discrete_x(j+1),discrete_y(j+1),theta(j),0); 
                lamda=results(j:j+1)';
                V_x=V_x+lamda*[t_a;t_b];
                V_y=V_y+lamda*[n_a;n_b];
            end
  
            %下一时刻例子的位置
            X(i,k+1)=X(i,k)+V_x*dt; 
            Y(i,k+1)=Y(i,k)+V_y*dt;
        end 
        k=k+1;
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