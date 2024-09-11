function [lamda,B] = solver(V,ctrl_x,ctrl_y,surf_x,surf_y,n ,length,theta)
%   涡元法求解
%   此处显示详细说明
syms X Y;
% A=zeros(2*n-1,2*n-2);
% RHS=zeros(2*n-1,1);
% use_linear=1;
% for i =1:2*n-2
%     RHS(i)=-dot(V,[cos(theta(i)),sin(theta(i))]);
%     for j=1:2*n-2
%         if j==i
%             A(i,j)=-1/2;
%         else
%             [u,w]=coefficient(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i)); 
%             A(i,j)=dot([u,w],[cos(theta(i)),sin(theta(i))]);
%         end
%     end
% end
% %库塔条件
% % A(k,:)=zeros(1,2*n-2);
% % RHS(k)=0;
% % A(k,n-1)=1;A(k,n)=1; 
% A(2*n-1,n-1)=1;A(2*n-1,n)=1; 
% %求解环量
% C=((A')*A)\((A')*RHS);
% % C=A\RHS;
% lamda=C(1:2*n-2);


A=zeros(2*n-1,2*n-1);
B=zeros(2*n-2,2*n-1);
RHS=zeros(2*n-1,1);
for i =1:2*n-2
    RHS(i)=-dot(V,[cos(theta(i)+pi/2),sin(theta(i)+pi/2)]);
    for j=1:2*n-2
   
        [n_a,n_b,t_a,t_b]=linear_cof(ctrl_x(i),ctrl_y(i),surf_x(j),surf_y(j),surf_x(j+1),surf_y(j+1),theta(j),theta(i),1); 
        A(i,j)=A(i,j)+n_a;
        A(i,j+1)=A(i,j+1)+n_b;            
        B(i,j)=B(i,j)+t_a;
        B(i,j+1)=B(i,j+1)+t_b;
    end
end
A(2*n-1,n-1)=1;A(2*n-1,n+1)=1; A(2*n-1,n)=2; 
lamda=A\RHS;

end

