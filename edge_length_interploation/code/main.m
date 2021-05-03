addpath('splsolver')

name1='./input_obj/bar_s.obj';
[points_old1,faces,~,~]=readObj(name1);
name2='./input_obj/bar_t.obj';
[points_old2,faces,~,~]=readObj(name2);

l1=edge_length(faces,points_old1);
l2=edge_length(faces,points_old2);
face_number=size(faces,1);
point_number=size(points_old1,1);
%A1,A2为原网格和目标网格的二面角
A1=compute_A(faces,points_old1);
A2=compute_A(faces,points_old2);
%point_to_faces为每个点相邻的面
%point_to_points为每个点相邻的点
[point_to_faces,point_to_points]=compute_point_neighbour(faces);

for iter=0:1:100
    t=iter/100;
    
    %插值二面角作为初始值
    A=(1-t)*A1+t*A2;
    %a为A的向量形式，方便计算
    %a(A_index)=A;
    [a,A_index]=A_to_a(faces,A);
    
    length=(1-t)*l1+t*l2;
    
    E1=energy_c(faces,point_number,a,point_to_faces,length,point_to_points,A_index)

    iteration=0;
    E_record=[];
    max_iteration=100;
    m=10;
    sk=zeros(max_iteration,size(a,1));
    yk=zeros(max_iteration,size(a,1));
    rho=zeros(size(a,1),1);
    
    d=zeros(3*face_number/2,1);
    g=zeros(3*face_number/2,1);
    
    %优化能量E,得到局部最优的二面角A
    %迭代所使用的是LBFGS方法
    while E1>1e-10 &&iteration<max_iteration
        iteration=iteration+1;
        if iteration==1
            %Mk为论文中的Mk
            Mk=compute_Mk_c(faces,point_number,A,point_to_faces,length,point_to_points);
            %gradient是二面角的梯度
            %g为gradient的向量形式
            gradient=compute_gradient(Mk,faces,A);
            g=gradient_to_g(faces,gradient,A_index);
            d=g;
        end
        
        E=E1;
        step=1;
        for i=1:150

            E1=energy_c(faces,point_number,a-d*step,point_to_faces,length,point_to_points,A_index);
            if E1<E
                a_temp=a;
                g_temp=g;
                a=a-d*step;
                A=atoA(faces,A_index,a);
                
                Mk=compute_Mk_c(faces,point_number,A,point_to_faces,length,point_to_points);
                gradient=compute_gradient(Mk,faces,A);
                g=gradient_to_g(faces,gradient,A_index);              
                             
                
                s=a-a_temp;
                y=g-g_temp;

                sk(iteration,:)=s;
                yk(iteration,:)=y;

                [sk,yk,rho,d]=LBFGS(g,iteration,m,sk,yk,rho);
                E_record=[E_record E1];
                E1
                break;
            else
                step=step/2;
            end
        end
    end
    save(sprintf("output/bar_%f.mat",t),"a");
    %画出能量变化
    figure
    plot(1:size(E_record,2),log10(E_record));
    xlabel('bar')
    drawnow;
    
    %M为相邻两点的局部坐标系的变换矩阵
    M=compute_M(A,faces,point_to_points,length,point_to_faces);
    %G1为点1的坐标架在世界坐标系下的值
    %G则为每个点的坐标架在世界坐标系下的值
    G1=eye(3);
    G=compute_G(G1,M,faces);
    G=orthogonalize_G(G);
    %local_coordinate为每个点的相邻点在其局部坐标系下的坐标
    local_coordinate=compute_local_coordinate(faces,point_to_points,length,point_to_faces,A);
    
    %重建出最终结果
    points=reconstruction(faces,local_coordinate,G,[0 0 0]);
    save_obj(points,faces,sprintf("output/bar_%f.obj",t));
    l3=edge_length(faces,points);
    sum(abs(length(:)-l3(:)));
end
