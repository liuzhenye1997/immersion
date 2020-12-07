function immersion(max_iteration,source_name,target_name,frames,Output_path,show_energy)
    format long
    addpath('splsolver')
    [points_old1,faces,~,~]=readObj(source_name);
    l1=edge_length(faces,points_old1);
    [points_old2,faces,~,~]=readObj(target_name);

    l2=edge_length(faces,points_old2);

    point_number=size(points_old1,1);

    fix=[1 point_number+1 point_number*2+1];
    free=1:3*point_number;
    free(fix)=[];

    l_target=l2;
    time=[];
    E_record=[];
    t_record=[];
    Time=[];
    Alpha_record=[];
    ZD_record=[];

    T1=clock;
    threshold=1e-29;
    dt=1/frames;
    for t=0.0:dt:1
        t
        l_target=t*l2+(1-t)*l1;
        points_temp=LSCM_sparse(l_target,faces);
        points_temp(:,3)=points_temp(:,3)+10^-10*(rand(size(points_temp(:,3)))-0.5);
     

        H=eye(size(points_old1,1)*3);
      
        energy_prev=Energy(faces,points_temp,l_target)

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
        m=rand(1)*13;
        I=sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,10^-m);
        H=hessian_matrix_sparse(points_temp,faces,l_target)+I;
        H1=H(free,free);
        if iteration==1
            solver = splsolver(H1, 'ldlt'); 
        end
        grad=compute_grad(points_temp,faces,l_target);
%         solver.refactorize(H);
%         x1 = solver.solve(H*points_temp(:)-grad(:));
%         x1=H1\(H1*points_temp(free)'-grad(free)');
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
%         alpha=100;
        E=Energy(faces,points_temp,l_target)
%         for i=1:100
%             p=points_temp;
%             p(2:end,:)=p(2:end,:)-alpha*d(2:end,:);
%             if Energy(faces,p,l_target)<E
%                 
%                 break;
%             else
%                 alpha=alpha/2;
%             end
%         end
        alpha=best_step(points_temp,faces,d,l_target);

        points_prev=points_temp;
        points_temp(free) = points_temp(free)-alpha*d(free);
        points_temp=reshape(points_temp,size(points_old1));
%         s=points_temp-points_prev;
%         y=compute_grad(points_temp,faces,l_target)-grad;
%         s=s(:);
%         y=y(:);
%         Hy=H*y;
%         H=H+(1+(y'*Hy)/(s'*y))*(s*s')/(s'*y)-((Hy)*s'+s*(Hy'))/(s'*y);
        
        e=Energy(faces,points_temp,l_target);
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
        
%         e=Energy(faces,points_temp,l_target);
        if e>2*E 
            e
            points_temp=reshape(points_prev,size(points_old1));
            break;
        end
        if alpha==0
            continue
        end
        if e<1e-29
            break
        end
        if mod(iteration,max_iteration)==0
%             figure 
%             plot(time1,log10(e_record));drawnow;
%             figure 
%             plot(time1,alpha_record);drawnow;
%             figure 
%             plot(time1,zd_record);drawnow;
            break;
        end
        end
        close
        if show_energy
            figure 
            plot(time1,log10(e_record));drawnow;
        end
%         figure 
%         plot(time1,alpha_record);drawnow;
%         figure 
%         plot(time1,zd_record);drawnow;
        objname=Output_path+sprintf("_%f.obj",t);
        matname=Output_path+sprintf("_%f.mat",t);
        save(matname,'points_temp');
        save_obj(points_temp,faces,objname);
        figure 
        trimesh(faces, points_temp(:,1), points_temp(:,2), points_temp(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
        drawnow; pause(0.001);
        l_temp=edge_length(faces,points_temp);
        energy=sum((l_temp(:).^2-l_target(:).^2).^2)
    end

    l_temp=edge_length(faces,points_temp);
    energy=sum((l_temp(:).^2-l2(:).^2).^2);
end

function grad=compute_grad(points_temp,faces,l_target)
%     grad=zeros(size(points_temp));
    l_temp=edge_length(faces,points_temp);
%     for i=1:size(faces,1)
%         for j=1:3
%             for k=1:3
%                 grad(faces(i,j),k)=grad(faces(i,j),k)+2*(l_temp(i,j)^2-l_target(i,j)^2)*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
%                 grad(faces(i,mod(j,3)+1),k)=grad(faces(i,mod(j,3)+1),k)-2*(l_temp(i,j)^2-l_target(i,j)^2)*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
%             end
%         end
%     end
    grad=compute_grad_c( points_temp, faces, l_target, l_temp);
end

function E=Energy(faces,points,l_target)
    l_temp=edge_length(faces,points);
    E=sum((l_temp(:).^2-l_target(:).^2).^2);
end
