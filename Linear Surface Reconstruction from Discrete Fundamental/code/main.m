addpath('splsolver')

name1='./input_obj/cat_s.obj';
[points_old1,faces,~,~]=readObj(name1);
name2='./input_obj/cat_t.obj';
[points_old2,faces,~,~]=readObj(name2);
l1=edge_length(faces,points_old1);
l2=edge_length(faces,points_old2);
face_number=size(faces,1);
point_number=size(points_old1,1);

%fr1为原网格的每个面上的坐标系框架,R1则为相邻面坐标系的变换矩阵,a1_1和a1_2则为三角形的三条边在
%fr1这个坐标系上的系数。q_fr1为将fr1从矩阵转换为四元数后的结果。
%对于目标网格，各变量的命名规则也相似。
fr1=compute_frame(points_old1,faces);  
R1=compute_R(faces,fr1);
[a1_1,a1_2]=compute_a(points_old1,faces,fr1);
q_fr1=zeros(face_number,4);
for i=1:face_number
    q_fr1(i,:)=Quaternion.matrix_to_q(fr1(3*i-2:3*i,:));
end

fr2=compute_frame(points_old2,faces);  
R2=compute_R(faces,fr2);
[a2_1,a2_2]=compute_a(points_old2,faces,fr2);
q_fr2=zeros(face_number,4);
for i=1:face_number
    q_fr2(i,:)=Quaternion.matrix_to_q(fr2(3*i-2:3*i,:));
end

%Q为将原网格个平面上的坐标系转换为目标网格的坐标系的变换矩阵对应的四元数
%Q=[cos(q_theta) N*sin(q_theta)]
Q=Quaternion.batch_multiplication(q_fr2,[q_fr1(:,1) -q_fr1(:,2) -q_fr1(:,3) -q_fr1(:,4)]);
q_theta=acos(Q(:,1));
N=Q(:,2:4)./sin(q_theta);

%R3为将R1转换为R2所用的变换矩阵
%下面的步骤是用于计算插值R1到R2的过程的旋转角和旋转轴
S=zeros(face_number,3,3,3);
q_R=zeros(face_number,3,4);
q_R1=zeros(face_number,3,4);
q_R2=zeros(face_number,3,4);
R3=zeros(face_number,3,3,3);
for i=1:face_number
    for j=1:3
        q_R1(i,j,:)=Quaternion.matrix_to_q(reshape(R1(i,j,:,:),3,3));
        q_R2(i,j,:)=Quaternion.matrix_to_q(reshape(R2(i,j,:,:),3,3));
        
        R3(i,j,:,:)=reshape(R2(i,j,:,:),3,3)*reshape(R1(i,j,:,:),3,3)';
        [U,D,V]=svd(reshape(R3(i,j,:,:),3,3));
        S(i,j,:,:)=V*D*V';
        R=U*V';
        q_R(i,j,:)=Quaternion.matrix_to_q(R);
    end
end
%q_R=[cos(q_R_theta) N_R*sin(q_R_theta)]
q_R_theta=acos(q_R(:,:,1));
N_R=zeros(face_number,3,3);
for i=1:face_number
    for j=1:3
        if(q_R_theta(i,j)~=0)
            N_R(i,j,:)=q_R(i,j,2:4)./sin(q_R_theta(i,j));
        else
            N_R(i,j,:)=[1 0 0];
        end
    end
end


norm_a1=sqrt(a1_1.^2+a1_2.^2);
norm_a2=sqrt(a2_1.^2+a2_2.^2);
N_a1_1=a1_1./norm_a1;
N_a1_2=a1_2./norm_a1;
N_a2_1=a2_1./norm_a2;
N_a2_2=a2_2./norm_a2;
a_trans=[N_a2_1.*N_a1_1+N_a1_2.*N_a2_2 N_a2_2.*N_a1_1-N_a1_2.*N_a2_1];
a_trans_theta=atan2(a_trans(:,4:6),a_trans(:,1:3));



for t=0:0.01:1
    %show_frame(points_old1,faces,fr);
    
    %计算稀疏矩阵。其中A为论文中能量Ef的系数矩阵，C为能量Ex中x，即点的坐标的系数矩阵，D为fr，即每个面上
    %的坐标系对应的系数矩阵
    q_Rt=zeros(size(q_R1));
    q_Rt(:,:,1)=cos(q_R_theta*t);
    q_Rt(:,:,2:4)=N_R.*sin(q_R_theta*t);
    R=zeros(size(R1));
   
    for j=1:3
         q=Quaternion.batch_multiplication(q_Rt(:,j,:),q_R1(:,j,:));
         m=Quaternion.batch_q_to_matrix(q);
         for i=1:face_number         
            R(i,j,:,:)=m(3*i-2:3*i,:);
        end
    end

    q_a=[cos(a_trans_theta*t) sin(a_trans_theta*t)];
    a1=((1-t)+t*norm_a2./norm_a1).*(q_a(:,1:3).*a1_1-q_a(:,4:6).*a1_2);
    a2=((1-t)+t*norm_a2./norm_a1).*(q_a(:,4:6).*a1_1+q_a(:,1:3).*a1_2);
    
    A=conpute_Ef_matrix(faces,R);
    [C,D]=conpute_Ex_matrix(point_number,faces,[a1,a2]);
    C=C';
    D=D';
    A=A';
    %w为权重。M为将[x,fr]当做一个变量后对应的系数矩阵，以求解能量的最小值
    w=1000;
    M=[C'*C -C'*D;-D'*C D'*D+w*(A'*A)];
    
    %插值计算某个面的坐标系，用来当做初始值
    Q4=[cos(q_theta*t) N.*sin(q_theta*t)];
    Q4=Quaternion.batch_multiplication(Q4,q_fr1);   
    fr=Quaternion.batch_q_to_matrix(Q4);
   
  
    for i=1:face_number
        fr(3*i-2:3*i,:)=fr(3*i-2:3*i,:)';
    end
    for i=face_number:-1:1
        fr(3*i-2:3*i,:)=fr(3*i-2:3*i,:)*fr(1:3,:)';
    end
    %固定一个面的三个点和这个面上坐标系，然后求解能量的最小值
    %xf=[x f],即为最终结果
    d=3;
    index_static=faces(1,:);
    index_var=1:point_number;
    index_var(index_static)=[];

    index_left=[index_var point_number+d+1:point_number+3*face_number];
    M_left=M(index_left,index_left);
    M_right=[M(index_var,index_static) M(index_var,point_number+1:point_number+d);M(point_number+d+1:point_number+3*face_number,index_static) M(point_number+d+1:point_number+3*face_number,point_number+1:point_number+d)];

    points=(1-t)*points_old1+t*points_old2;
    l=(1-t)*l1+t*l2;
    theta=acos((l(1,1)^2+l(1,3)^2-l(1,2)^2)/(2*l(1,1)*l(1,3)));
    p=[0 0 0;l(1,1) 0 0;l(1,3)*cos(theta) l(1,3)*sin(theta) 0];
    xf1=[p;eye(3)];
    xf2=-M_left\(M_right*xf1);
    xf=[points;fr];
    xf(index_left,:)=xf2;
    xf(index_static,:)=p;
    

    objname=sprintf("output\\cat_%f.obj",t);
    points_temp=xf(1:point_number,:);
    save_obj(points_temp,faces,objname);
    l_temp=edge_length(faces,points_temp);
    
    sum(abs(l_temp(:)-l(:)))
end

