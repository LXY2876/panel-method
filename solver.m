function [lamda,B,D,E] = solver(V,ctrl_x,ctrl_y,surf_x,surf_y,n ,theta,Flag)
%   涡元法求解
%   此处显示详细说明
syms X Y;
Vinj_V=1;Vsuc_V=0.5;
 
 
N_points=length(Flag);%控制面个数
N_gama=sum(abs(Flag))+1+4;%未知量个数
A=zeros(N_gama,N_gama);
A2=zeros(N_gama,N_gama);
B=zeros(N_points,N_gama);%求解每个控制点的速度的矩阵
E=zeros(N_points,N_gama);%校验的矩阵
D=zeros(N_points,N_gama);%求解每个控制点的值的矩阵
RHS=zeros(N_gama,1);
kuta_index=[];
index_row=1;
for i=1:N_points

    if Flag(i)==1
        
        RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
        index_col=1;
        for j=1:N_points  
            if Flag(j)==1
                [n_a,n_b,t_a,t_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;            
                B(i,index_col)=B(i,index_col)+t_a;
                B(i,index_col+1)=B(i,index_col+1)+t_b;
                E(i,index_col)=E(i,index_col)+n_a;
                E(i,index_col+1)=E(i,index_col+1)+n_b;
            else 
                [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A2(index_row,index_col)=A2(index_row,index_col)+n_a;
                A2(index_row,index_col+1)=A2(index_row,index_col+1)+n_b;   
                A2(index_row,index_col+2)=A2(index_row,index_col+2)+n_c;   
                B(i,index_col)=B(i,index_col)+t_a;
                B(i,index_col+1)=B(i,index_col+1)+t_b;
                B(i,index_col+2)=B(i,index_col+2)+t_c;
                E(i,index_col)=E(i,index_col)+n_a;
                E(i,index_col+1)=E(i,index_col+1)+n_b;
                E(i,index_col+2)=E(i,index_col+2)+n_c;
            end
            if j==i
                D(i,index_col)=0.5;
                D(i,index_col+1)=0.5;
            end
            if Flag(j)~=Flag(mod(j,N_points)+1)
                index_col=index_col+1;
            end   

            index_col=index_col+abs(Flag(j));
        end
    elseif abs(Flag(i))==2
        if Flag(i)>0
            RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])+Vinj_V*norm(V);  %法向约束   
        else
            RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])-Vsuc_V*norm(V);
        end
        RHS(index_row+1)=-dot(V,[cos(theta(i)),sin(theta(i))]); %切向约束
        index_col=1;
        for j=1:N_points  
            if Flag(j)==1
                [n_a,n_b,t_a,t_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;            
                A(index_row+1,index_col)=A(index_row+1,index_col)+t_a;
                A(index_row+1,index_col+1)=A(index_row+1,index_col+1)+t_b;            
                B(i,index_col)=B(i,index_col)+n_a;
                B(i,index_col+1)=B(i,index_col+1)+n_b;
                E(i,index_col)=E(i,index_col)+t_a;
                E(i,index_col+1)=E(i,index_col+1)+t_b;
            else 
                [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A2(index_row,index_col)=A2(index_row,index_col)+n_a;
                A2(index_row,index_col+1)=A2(index_row,index_col+1)+n_b;   
                A2(index_row,index_col+2)=A2(index_row,index_col+2)+n_c;   
                A2(index_row+1,index_col)=A2(index_row+1,index_col)+t_a;
                A2(index_row+1,index_col+1)=A2(index_row+1,index_col+1)+t_b;   
                A2(index_row+1,index_col+2)=A2(index_row+1,index_col+2)+t_c;  
                B(i,index_col)=B(i,index_col)+n_a;
                B(i,index_col+1)=B(i,index_col+1)+n_b;
                B(i,index_col+2)=B(i,index_col+2)+n_c;
                E(i,index_col)=E(i,index_col)+t_a;
                E(i,index_col+1)=E(i,index_col+1)+t_b;
                E(i,index_col+2)=E(i,index_col+2)+t_c;
            end
            if j==i
                 
                D(i,index_col+1)=1;
            end
            if Flag(j)~=Flag(mod(j,N_points)+1)
                index_col=index_col+1;
            end

            index_col=index_col+abs(Flag(j));

        end
    end

    index_row=index_row+abs(Flag(i));
    if Flag(i)~=Flag(mod(i,N_points)+1)
        
        kuta_index=[kuta_index;i,index_row];
        index_row=index_row+1;
    end
end
% kuta_index=[kuta_index;N_gama-n,N_gama];
%尾缘库塔条件
RHS(N_gama)=-dot(V,[cos(theta(N_points-n+1)),sin(theta(N_points-n+1))])-dot(V,[cos(theta(N_points-n+2)),sin(theta(N_points-n+2))]);
index_row=N_gama;
index_col=1;
for j=1:N_points  
    
    if Flag(j)==1
        [~,~,t_a,t_b]=linear_cof(ctrl_x(N_points-n+1),ctrl_y(N_points-n+1),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(N_points-n+1)); 
        A(index_row,index_col)=A(index_row,index_col)+t_a;
        A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;            
        [~,~,t_a,t_b]=linear_cof(ctrl_x(N_points-n+2),ctrl_y(N_points-n+2),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(N_points-n+2)); 
        A(index_row,index_col)=A(index_row,index_col)+t_a;
        A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;           
    else 
        [~,~,~,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(N_points-n+1),ctrl_y(N_points-n+1),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(N_points-n+1)); 
        A(index_row,index_col)=A(index_row,index_col)+t_a;
        A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;   
        A(index_row,index_col+2)=A(index_row,index_col+2)+t_c;   
        [~,~,~,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(N_points-n+2),ctrl_y(N_points-n+2),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(N_points-n+2)); 
        A(index_row,index_col)=A(index_row,index_col)+t_a;
        A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;   
        A(index_row,index_col+2)=A(index_row,index_col+2)+t_c;   
    end
    if Flag(j)~=Flag(mod(j,N_points)+1)
        index_col=index_col+1;
    end            
    index_col=index_col+abs(Flag(j));
end

% A(N_gama,N_gama-n+2)=1;A(N_gama,N_gama-n)=1; A(N_gama,N_gama-n+1)=2; 
A=A+A2;
%交界处切向速度连续
for k= 1:4
    index_col=1;
    i=kuta_index(k,1);
    index_row=kuta_index(k,2);
    if Flag(i)*Flag(i+1)<0
        dircetion=-1;
    else
        dircetion=1;
    end
    if abs(Flag(i))==2  
%         RHS(index_row)=dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])-dircetion*dot(V,[cos(theta(i+1)),sin(theta(i+1))]);
        jet_ind=i;
        wall_ind=i+1;    
    else        
        jet_ind=i+1;
        wall_ind=i;
    end
    RHS(index_row)=-dot(V,[cos(theta(wall_ind)),sin(theta(wall_ind))])+dircetion*dot(V,[cos(theta(jet_ind)+pi/2),sin(theta(jet_ind)+pi/2)]);
    for j=1:N_points  
        if Flag(j)==1
            [~,~,t_a,t_b]=linear_cof(ctrl_x(wall_ind),ctrl_y(wall_ind),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(wall_ind)); 
            A(index_row,index_col)=A(index_row,index_col)+t_a;
            A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;            
            [n_a,n_b,~,~]=linear_cof(ctrl_x(jet_ind),ctrl_y(jet_ind),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(jet_ind)); 
            A(index_row,index_col)=A(index_row,index_col)-dircetion*n_a;
            A(index_row,index_col+1)=A(index_row,index_col+1)-dircetion*n_b;           
        else 
            [~,~,~,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(wall_ind),ctrl_y(wall_ind),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(wall_ind)); 
            A(index_row,index_col)=A(index_row,index_col)+t_a;
            A(index_row,index_col+1)=A(index_row,index_col+1)+t_b;   
            A(index_row,index_col+2)=A(index_row,index_col+2)+t_c;   
            [n_a,n_b,n_c,~,~,~]=Quadratic_source_cof(ctrl_x(jet_ind),ctrl_y(jet_ind),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(jet_ind)); 
            A(index_row,index_col)=A(index_row,index_col)-dircetion*n_a;
            A(index_row,index_col+1)=A(index_row,index_col+1)-dircetion*n_b;   
            A(index_row,index_col+2)=A(index_row,index_col+2)-dircetion*n_c;   
        end
        if Flag(j)~=Flag(mod(j,N_points)+1)
            index_col=index_col+1;
        end            
        index_col=index_col+abs(Flag(j));
    end
end


% lamda=A\RHS;
[P,R,C] = equilibrate(A);
B1=R*P*A*C;
d=R*P*RHS;
y=B1\d;
lamda=C*y;
end
 