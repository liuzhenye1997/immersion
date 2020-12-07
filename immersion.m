function immersion(max_iteration,source_name,target_name,frames,Output_path,show_energy,show_result)
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

    threshold=1e-19;
    dt=1/frames;
    for t=0.0:dt:1
        t
        l_target=t*l2+(1-t)*l1;
        points_temp=LSCM_sparse(l_target,faces);
        points_temp(:,3)=points_temp(:,3)+10^-10*(rand(size(points_temp(:,3)))-0.5);
     
        energy_prev=Energy(faces,points_temp,l_target)

        e_record=[];    
 
  
        t1=clock;
        time1=0;
        iteration=0;

        while 1
            iteration=iteration+1;
        if iteration==1
            [row,col,~,H]=svd_positive_hessian_matrix_sparse_c(points_temp,faces,l_target);
            H=H+sparse(1:size(points_old1,1)*3,1:size(points_old1,1)*3,10^-15);
            H_index=ij2nzIdxs(H, uint64(row(:)), uint64(col(:)));
            [M,N]=find(H);
            k=unique([find(sum(M(:)==fix,2)) ,find(sum(N(:)==fix,2)) ]);
            H1_index=1:size(M,1);
            H1_index(k)=[];
            H1=H(free,free);
        else
            value=svd_positive_hessian_matrix_value_c(points_temp,faces,l_target);
            H_val=myaccumarray(H_index, value);
            update_matrix(H1,H_val,H1_index)

        end
            if iteration==1
                solver = splsolver(H1, 'llt'); 
            end
            grad=compute_grad(points_temp,faces,l_target);

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
            alpha=best_step(points_temp,faces,d,l_target);

            points_temp(free) = points_temp(free)-alpha*d(free);
            points_temp=reshape(points_temp,size(points_old1));


            e=Energy(faces,points_temp,l_target);
            time1=[time1 etime(clock,t1)];
     
            if mod(iteration,max_iteration)==0
                break;
            end

            if alpha==0
                continue
            end
            if e<1e-29
                break
            end
          
        end

        if show_energy
            figure 
            plot(time1,log10(e_record));drawnow;
        end

        objname=Output_path+sprintf("_%f.obj",t);
        matname=Output_path+sprintf("_%f.mat",t);
        save(matname,'points_temp');
        save_obj(points_temp,faces,objname);
        if show_result
            figure 
            trimesh(faces, points_temp(:,1), points_temp(:,2), points_temp(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
            drawnow; pause(0.001);
        end
        l_temp=edge_length(faces,points_temp);
        energy=sum((l_temp(:).^2-l_target(:).^2).^2)
    end
end

function grad=compute_grad(points_temp,faces,l_target)
    l_temp=edge_length(faces,points_temp);
    grad=compute_grad_c( points_temp, faces, l_target, l_temp);
end

function E=Energy(faces,points,l_target)
    l_temp=edge_length(faces,points);
    E=sum((l_temp(:).^2-l_target(:).^2).^2);
end
