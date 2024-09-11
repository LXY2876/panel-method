function [lamda,B] = solver(V,ctrl_x,ctrl_y,surf_x,surf_y,n ,theta,Flag)
%   涡元法求解
%   此处显示详细说明
syms X Y;
Vinj_V=0.1;Vsuc_V=0.05;
 
 
N_points=length(Flag);%控制点个数
N_gama=sum(abs(Flag))+1;%未知量个数
A=zeros(N_gama,N_gama);
A2=zeros(N_gama,N_gama);
B=zeros(N_gama,N_gama);
RHS=zeros(N_gama,1);
 
index_row=1;
for i=1:N_points

    if Flag(i)==1
        RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
        index_col=1;
        for j=1:N_points  
            if Flag(j)==1
                [u_a,w_a,u_b,w_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+dot([u_a,w_a],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
                A(index_row,index_col+1)=A(index_row,index_col+1)+dot([u_b,w_b],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);            
            else
                [u_a,w_a,u_b,w_b,u_c,w_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A2(index_row,index_col)=A2(index_row,index_col)+dot([u_a,w_a],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
                A2(index_row,index_col+1)=A2(index_row,index_col+1)+dot([u_b,w_b],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);   
                A2(index_row,index_col+2)=A2(index_row,index_col+2)+dot([u_c,w_c],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);   
            end
            index_col=index_col+abs(Flag(j));
        end
    elseif abs(Flag(i))==2
        if Flag(i)>0
            RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])+Vinj_V*norm(V);  %法向约束   
        else
            RHS(index_row)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])-Vsuc_V*norm(V);
        end
        RHS(index_row+1)=-dot(V,[cos(theta(i)),sin(theta(i))]);
        index_col=1;
        for j=1:N_points  
            if Flag(j)==1
                [u_a,w_a,u_b,w_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+dot([u_a,w_a],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
                A(index_row,index_col+1)=A(index_row,index_col+1)+dot([u_b,w_b],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);            
                A(index_row+1,index_col)=A(index_row+1,index_col)+dot([u_a,w_a],[cos(theta(i)),sin(theta(i))]);
                A(index_row+1,index_col+1)=A(index_row+1,index_col+1)+dot([u_b,w_b],[cos(theta(i)),sin(theta(i))]);            
            
            elseif abs(Flag(i))==2
                [u_a,w_a,u_b,w_b,u_c,w_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A2(index_row,index_col)=A2(index_row,index_col)+dot([u_a,w_a],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
                A2(index_row,index_col+1)=A2(index_row,index_col+1)+dot([u_b,w_b],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);   
                A2(index_row,index_col+2)=A2(index_row,index_col+2)+dot([u_c,w_c],[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);   
                A2(index_row+1,index_col)=A2(index_row+1,index_col)+dot([u_a,w_a],[cos(theta(i)),sin(theta(i))]);
                A2(index_row+1,index_col+1)=A2(index_row+1,index_col+1)+dot([u_b,w_b],[cos(theta(i)),sin(theta(i))]);   
                A2(index_row+1,index_col+2)=A2(index_row+1,index_col+2)+dot([u_c,w_c],[cos(theta(i)),sin(theta(i))]);  
            else
                continue;
            end
            index_col=index_col+abs(Flag(j));
        end
    else

    
    end
    index_row=index_row+abs(Flag(i));
end
 
%尾缘库塔条件
A(N_gama,N_gama-n+1)=1;A(N_gama,N_gama-n-1)=1; A(N_gama,N_gama-n)=2; 
 


% lamda=A\RHS;
[P,R,C] = equilibrate(A+A2);
B=R*P*A*C;
d=R*P*RHS;
y=B\d;
lamda=C*y;
end
 