function [lamda,B,D,E] = solver(V,ctrl_x,ctrl_y,surf_x,surf_y,theta,Flag)
%   涡元法求解
%   此处显示详细说明
syms X Y;
Vinj_V=1;Vsuc_V=1; 

N_panels=length(Flag);
N_gama=sum(abs(Flag))+4;
A=zeros(N_gama,N_gama);
B=zeros(N_panels,N_gama);
D=zeros(N_panels,N_gama);
E=zeros(N_panels,N_gama);%校验的矩阵
RHS=zeros(N_gama,1);
index_row=1;
kuta_index=[];
for i =1:N_panels
     
    index_col=1;
    if abs(Flag(i))==1
        for j=1:N_panels
            if abs(Flag(j))==1
                [n_a,n_b,t_a,t_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;        
                B(i,index_col)=B(i,index_col)+t_a;
                B(i,index_col+1)=B(i,index_col+1)+t_b;
                E(i,index_col)=E(i,index_col)+n_a;
                E(i,index_col+1)=E(i,index_col+1)+n_b;
            elseif abs(Flag(j))==2
                [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;    
                A(index_row,index_col+2)=A(index_row,index_col+2)+n_c;   
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
            if Flag(j)~=Flag(mod(j,N_panels)+1)
                index_col=index_col+1;
            end   
            index_col=index_col+abs(Flag(j));
        end
    else
        if Flag(i)>0
            RHS(index_row)=Vinj_V*norm(V);  %法向约束   
        else
            RHS(index_row)=-Vsuc_V*norm(V);
        end
        for j=1:N_panels
            if abs(Flag(j))==1
                [n_a,n_b,t_a,t_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;   
                A(index_row+1,index_col)=A(index_row+1,index_col)+t_a;
                A(index_row+1,index_col+1)=A(index_row+1,index_col+1)+t_b;   
                B(i,index_col)=B(i,index_col)+n_a;
                B(i,index_col+1)=B(i,index_col+1)+n_b;
                E(i,index_col)=E(i,index_col)+t_a;
                E(i,index_col+1)=E(i,index_col+1)+t_b;
            elseif abs(Flag(j))==2
                [n_a,n_b,n_c,t_a,t_b,t_c]=Quadratic_source_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
                A(index_row,index_col)=A(index_row,index_col)+n_a;
                A(index_row,index_col+1)=A(index_row,index_col+1)+n_b;    
                A(index_row,index_col+2)=A(index_row,index_col+2)+n_c;    
                A(index_row+1,index_col)=A(index_row+1,index_col)+t_a;
                A(index_row+1,index_col+1)=A(index_row+1,index_col+1)+t_b;    
                A(index_row+1,index_col+2)=A(index_row+1,index_col+2)+t_c;   
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
            if Flag(j)~=Flag(mod(j,N_panels)+1)              
                index_col=index_col+1; 
            end   
            index_col=index_col+abs(Flag(j));
        end
    end
    index_row=index_row+abs(Flag(i));
    if Flag(i)~=Flag(mod(i,N_panels)+1)
        
        kuta_index=[kuta_index;i,i+1,index_row];
        index_row=index_row+1;
    end
end
% A(2*n-1,n-1)=1;A(2*n-1,n+1)=1; A(2*n-1,n)=2; 
kuta_index(4,2)=1;
%交界处切向速度连续
for k= 1:4
    index_col=1;
  
    index_row=kuta_index(k,3);
    if Flag(kuta_index(k,1))*Flag(kuta_index(k,2))<0
        dircetion=-1;
    else
        dircetion=1;
    end
    if abs(Flag(kuta_index(k,1)))==2  
%         RHS(index_row)=dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)])-dircetion*dot(V,[cos(theta(i+1)),sin(theta(i+1))]);
        jet_ind=kuta_index(k,1);
        wall_ind=kuta_index(k,2);    
    else        
        jet_ind=kuta_index(k,2);  
        wall_ind=kuta_index(k,1);
    end

    for j=1:N_panels  
        if abs(Flag(j))==1
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
        if Flag(j)~=Flag(mod(j,N_panels)+1)
            index_col=index_col+1;
        end            
        index_col=index_col+abs(Flag(j));
    end
end


[P,R,C] = equilibrate(A);
B1=R*P*A*C;
d=R*P*RHS;
y=B1\d;
lamda=C*y;


end

