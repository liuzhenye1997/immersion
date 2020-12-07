function immersion3D(max_iteration,target_name,Output_path,show_energy,show_result,energy_type,initialization_obj)
    format long
    addpath('splsolver')
    if energy_type~=1 && energy_type~=2 && energy_type~=3 && energy_type~=4 && energy_type~=5
        error("energy_type应为1,2,3,4或者5");
    end
    [points_old1,faces,~,~]=readObj(target_name);

    l1=edge_length(faces,points_old1);

    point_number=size(points_old1,1);

    fix=[1 point_number+1 point_number*2+1];
    free=1:3*point_number;
    free(fix)=[];

    l_target=l1;
    time=[];
    E_record=[];
    t_record=[];
    Time=[];
    Alpha_record=[];
    ZD_record=[];

    T1=clock;
    threshold=1e-19;
    for t=0
        l_target=l1;
        if nargin==6
            points_temp=shape_from_metric(10,l_target,faces);
        elseif nargin==7
            [points_temp,~]=readObj(char(initialization_obj));
        end
    
        H=eye(size(points_old1,1)*3);
        I=sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,1e-6);
        energy_prev=Energy(faces,points_temp,l_target,energy_type)

        E_record=[E_record energy_prev];
        e_record=[];    
        e_record=[e_record energy_prev];
        alpha_record=[];
        alpha_record=[alpha_record 0];
        Alpha_record=[Alpha_record 0];
        ZD_record=[ZD_record 0];
        zd_record=0;
        time=[time 0];
        t_record=[t_record t];
        t1=clock;
        Time=[Time etime(clock,T1)];
        time1=0;
        iteration=0;

        while 1
            iteration=iteration+1;
        if iteration==1
            if energy_type==1
                [row,col,~,H]=svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                H=H+sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,10^-15);
                H_index=ij2nzIdxs(H, uint64(row(:)), uint64(col(:)));
                [M,N]=find(H);
                k=unique([find(sum(M(:)==fix,2)) ,find(sum(N(:)==fix,2)) ]);
                H1_index=1:size(M,1);
                H1_index(k)=[];        
            elseif energy_type==2
                [~,~,~,H]=sqrt_energy_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                H=H+sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,10^-5);
                H1=H(free,free);
            elseif energy_type==3
                [~,~,~,H]=energy3_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                H1=H(free,free);
            elseif energy_type==4
                [~,~,~,H]=energy4_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                H1=H(free,free); 
            elseif energy_type==5
                [~,~,~,H]=energy3_1_svd2_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                 H1=H(free,free);
            end
        else
            if energy_type==1
                value=svd_positive_hessian_matrix_value_c(points_temp,faces,l_target);
                H_val=myaccumarray(H_index, value);
                update_matrix_single_thread(H1,H_val,H1_index)
            elseif energy_type==2
               [~,~,~,H]=sqrt_energy_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);              
               H=H++sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,10^-5);
               H1=H(free,free);
            elseif energy_type==3
                [~,~,~,H]=energy3_svd2_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
%                [~,~,~,H]=energy3_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
               H1=H(free,free);
            elseif energy_type==4
                [~,~,~,H]=energy4_svd2_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
%                [~,~,~,H]=energy4_svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
               H1=H(free,free); 
            elseif energy_type==5
                 [~,~,~,H]=energy3_1_svd2_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
                 H1=H(free,free);
            end
   
        end
            if iteration==1
                if energy_type~=2
                    solver = splsolver(H1, 'llt'); 
                else
                    solver = splsolver(H1, 'ldlt'); 
                end
            end
            grad=compute_grad(points_temp,faces,l_target,energy_type);

            x1 =solver.refactor_solve(H1,H1*points_temp(free)'-grad(free)');
            d=x1-points_temp(free)';
            z=grad(free);
            zd=z*d(:);
            if zd<0
                d=-d;
            end 
            if sum(abs(z(:)))/size(z,1)<threshold
                sum(abs(z(:)))/size(z,1)
                break
            end
            D=zeros(3*point_number,1);
            D(free)=d;
            d=reshape(D,size(points_temp));
            E=Energy(faces,points_temp,l_target,energy_type)
            if energy_type==1
                alpha=best_step(points_temp,faces,d,l_target);
            else
                alpha=100;          
                for i=1:100
                    p=points_temp;
                    p(2:end,:)=p(2:end,:)-alpha*d(2:end,:);
                    if Energy(faces,p,l_target,energy_type)<E
                        break;
                    else
                        alpha=alpha/2;
                    end
                end
            end

            points_prev=points_temp;
            points_temp(free) = points_temp(free)-alpha*d(free);
            points_temp=reshape(points_temp,size(points_old1));


            e=Energy(faces,points_temp,l_target,energy_type);
            E_record=[E_record e];
            e_record=[e_record e];
            time=[time etime(clock,t1)];
            time1=[time1 etime(clock,t1)];
            t_record=[t_record t];
            Time=[Time etime(clock,T1)];
            Alpha_record=[Alpha_record alpha];
            alpha_record=[alpha_record alpha];
            ZD_record=[ZD_record zd];
            zd_record=[zd_record zd];
            if mod(iteration,max_iteration)==0
                break;
            end
            if e>2*E 
                e
                points_temp=reshape(points_prev,size(points_old1));
                break;
            end
            if alpha==0
                continue
            end
            if e<1e-20
                break
            end
          
        end
%         close
        if show_energy
            figure 
            plot(time1,log10(e_record));drawnow;
        end
%         figure 
%         plot(time1,alpha_record);drawnow;
%         figure 
%         plot(time1,zd_record);drawnow;
        objname=Output_path+sprintf("_%d.obj",max_iteration);
        matname=Output_path+sprintf("_%d.mat",max_iteration);
        save(matname,'points_temp');
        save_obj(points_temp,faces,objname);
        if show_result
            figure 
            trimesh(faces, points_temp(:,1), points_temp(:,2), points_temp(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
            drawnow; pause(0.001);
        end
    end
  
    mesh = Triangular_mesh(faces,l1);
    [max_abs_rel_err,min_rel_err,max_rel_err,l2_rel_err,max_log_err,min_log_err]=evaluate(mesh.f_number,mesh.v_dst,mesh.v_src,points_temp,mesh.length);
    fprintf("最大相对误差为%f,L2误差为%f",max_abs_rel_err,l2_rel_err);
end



function grad=compute_grad(points_temp,faces,l_target,energy_type)
    l_temp=edge_length(faces,points_temp);
    if energy_type==1
        grad=compute_grad_c( points_temp, faces, l_target, l_temp);
    elseif energy_type==2
        grad=sqrt_energy_compute_grad(points_temp,faces,l_target,l_temp);
    elseif energy_type==3
        grad=energy3_compute_grad(points_temp,faces,l_target);  
    elseif energy_type==4
        grad=energy4_compute_grad(points_temp,faces,l_target);
    elseif energy_type==5
        weight=1./l_target.^2;
        grad=energy3_1_compute_grad(points_temp,faces,l_target,weight);
    end
end

function grad=sqrt_energy_compute_grad(points_temp,faces,l_target,l_temp)
    grad=zeros(size(points_temp));
    for i=1:size(faces,1)
        for j=1:3
            for k=1:3
                grad(faces(i,j),k)=grad(faces(i,j),k)+sign(l_temp(i,j)^2-l_target(i,j)^2)/sqrt(abs(l_temp(i,j)^2-l_target(i,j)^2))*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
                grad(faces(i,mod(j,3)+1),k)=grad(faces(i,mod(j,3)+1),k)-sign(l_temp(i,j)^2-l_target(i,j)^2)/sqrt(abs(l_temp(i,j)^2-l_target(i,j)^2))*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
            end
        end
    end
end

function E=Energy(faces,points,l_target,energy_type)
    l_temp=edge_length(faces,points);
    if energy_type==1
        E=sum((l_temp(:).^2-l_target(:).^2).^2);
    elseif energy_type==2
        E=sum(abs(l_temp(:).^2-l_target(:).^2).^(1/2));
    elseif energy_type==3
        E=sum((l_temp(:)-l_target(:)).^2);
    elseif energy_type==4
        E=sum(exp((l_temp(:)-l_target(:)).^2)-1);
    elseif energy_type==5
        weight=1./l_target.^2;
        E=sum(weight(:).*(l_temp(:)-l_target(:)).^2);
    end
end

function grad=energy4_compute_grad(points_temp,faces,l_target)
    grad=zeros(size(points_temp));
    l_temp=edge_length(faces,points_temp);
    for i=1:size(faces,1)
        for j=1:3
            for k=1:3
                g=2*exp((l_temp(i,j)-l_target(i,j))^2)*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k))*(l_temp(i,j)-l_target(i,j))/l_temp(i,j);
                grad(faces(i,j),k)=grad(faces(i,j),k)+g;
                grad(faces(i,mod(j,3)+1),k)=grad(faces(i,mod(j,3)+1),k)-g;
            end
        end
    end
end
