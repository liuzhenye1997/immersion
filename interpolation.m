function interpolation(name1,name2,energy_type,p)
    max_iteration=1000;
    [points_old1,faces,~,~]=readObj(name1);
    l1=edge_length(faces,points_old1);

    [points_old2,faces,~,~]=readObj(name2);
    static=faces(1,:);
    % points_old2(static,:)=points_old1(static,:);
    l2=edge_length(faces,points_old2);
    % energy=sum((l1(:).^2-l2(:).^2).^2)
    point_number=size(points_old1,1);
    points_temp=points_old1;%+(rand(size(points_old1))-0.5)*0.01;
    % points_temp(2:end,:)=points_temp(2:end,:)+10^-8*(rand(size(points_temp(2:end,:)))-0.5)*0.01;
    points_temp_backup=points_temp;
    l_temp=edge_length(faces,points_temp);
    energy=sum((l_temp(:).^2-l2(:).^2).^2);
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
    % w=0.001;
    constant=1e3
    w=(sum(l1(:).^2+l2(:).^2))/(size(faces,1)*3*2)/constant;
    T1=clock;
    threshold=1e-19;
    for t=0:0.01:1
        objname=sprintf("SQP\\result\\%d_%d\\%s_%f.obj",energy_type,p,name1(1:end-13),t);
    %         if exist(objname,"file")
    %             continue;
    %         end
        l_target=sqrt(t*l2.^2+(1-t)*l1.^2);
        %     points_temp=LSCM_sparse(l_target,faces);
        %     name=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\max_iteration=2000\\blue_monster_%f.obj",t))
        %         name=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\LSCM_拟牛顿_阈值3e-11\\blue_monster_%f.obj",t))
        %         name=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\3\\blend_I_Newton_LSCM_continue3\\I=1e-0——1e-13\\max_iteration=20000\\blue_monster_%f.obj",t))
        %     name=char(sprintf("result\\I=1e-0——1e-13\\max_iteration=100\\blue_monster_%f.obj",t));
        %     [points_temp,faces,~,~]=readObj(name);
        %     points_temp=LSCM_sparse(l_target,faces);
        %     if t==0
        %         points_temp=points_old1;
        name=char(sprintf("%s_%f.obj",name1(1:end-13),t))
        [points_temp,faces,~,~]=readObj(name);
        initial_points=points_temp;


        energy_prev=Energy(faces,points_temp,l_target,energy_type,p,w,initial_points)



        %     H=eye(size(points_old1,1)*3);



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
        if energy_prev<1e-27
            objname=sprintf("SQP\\result\\%d_%d\\%s_%f.obj",energy_type,p,name1(11:end-13),t)
            matname=sprintf("SQP\\result\\%d_%d\\%s_%f.mat",energy_type,p,name1(11:end-13),t);
            save(matname,'points_temp');
            save_obj(points_temp,faces,objname);
            continue
        end
        while 1
            iteration = iteration + 1;
            if iteration == 1
                if energy_type == 5 || energy_type == 4 || energy_type == 3 || energy_type == 2 || energy_type == 1
                    [row, col, ~, H] = svd_positive_hessian_matrix_c(points_temp, faces, l_target, energy_type);
                elseif energy_type == 6
                    [row, col, ~, H] = svd_positive_hessian_matrix_c(points_temp, faces, l_target, energy_type, p);
                 elseif energy_type == 7
                    [row, col, ~, H] = svd_positive_hessian_matrix_c(points_temp, faces, l_target, 1);
                end
                H = H + sparse(1:size(points_old1, 1)*3, 1:size(points_old1, 1)*3, 10^-15);
                H_index = ij2nzIdxs(H, uint64(row(:)), uint64(col(:)));
                [M, N] = find(H);
                k = unique([find(sum(M(:) == fix, 2)), find(sum(N(:) == fix, 2))]);
                H1_index = 1:size(M, 1);
                H1_index(k) = [];
                H1 = H(free, free);
                if energy_type==7
                    H1=H1+sparse(1:size(H1,1),1:size(H1,1),2*w);
                end
            else
                if energy_type == 5 || energy_type == 4 || energy_type == 3 || energy_type == 2 || energy_type == 1
                    value = svd_positive_hessian_matrix_value_c(points_temp, faces, l_target, energy_type);
                elseif energy_type == 6
                    value = svd_positive_hessian_matrix_value_c(points_temp, faces, l_target, energy_type, p);
                elseif energy_type == 7
                    value = svd_positive_hessian_matrix_value_c(points_temp, faces, l_target,1);
                end

                H_val = myaccumarray(H_index, value);
                update_matrix_single_thread(H1, H_val, H1_index)
                if energy_type==7
                    H1=H1+sparse(1:size(H1,1),1:size(H1,1),2*w);
                end
            end


            if iteration == 1
                if energy_type ~= 2
                    solver = splsolver(H1, 'ldlt');
                else
                    solver = splsolver(H1, 'ldlt');
                end
            end
            if energy_type < 6
                grad = compute_grad(points_temp, faces, l_target, energy_type);
            elseif energy_type==6
                grad = compute_grad(points_temp, faces, l_target, energy_type, p);
            elseif energy_type==7
                grad = compute_grad(points_temp, faces, l_target, energy_type, p,w,initial_points);
    %             grad1 = compute_grad(points_temp, faces, l_target, 1);
            end


            x1 = solver.refactor_solve(H1, H1*points_temp(free)'-grad(free)');
            d = x1 - points_temp(free)';
            z = grad(free);
            zd = z * d(:);
            if zd < 0
                d = -d;
            end
            if sum(abs(z(:))) / size(z, 1) < threshold
                sum(abs(z(:))) / size(z, 1)
                break
            end
            D = zeros(3*point_number, 1);
            D(free) = d;
            d = reshape(D, size(points_temp));
            E = Energy(faces, points_temp, l_target, energy_type, p,w,initial_points)
            if energy_type == 1
                alpha = best_step(points_temp, faces, d, l_target);
            else
                alpha = 10;
                for i = 1:100
                    points_temp1 = points_temp;
                    points_temp1(2:end, :) = points_temp1(2:end, :) - alpha * d(2:end, :);
                    if Energy(faces, points_temp1, l_target, energy_type,  p,w,initial_points) < E
                        break;
                    else
                        alpha = alpha / 2;
                    end
                end
            end

            points_prev = points_temp;
            points_temp(free) = points_temp(free) - alpha * d(free);
            points_temp = reshape(points_temp, size(points_old1));


            e = Energy(faces, points_temp, l_target, energy_type,  p,w,initial_points);
            E_record = [E_record, e];
            e_record = [e_record, e];
            time = [time, etime(clock, t1)];
            time1 = [time1, etime(clock, t1)];
            t_record = [t_record, t];
            Alpha_record = [Alpha_record, alpha];
            alpha_record = [alpha_record, alpha];
            ZD_record = [ZD_record, zd];
            zd_record = [zd_record, zd];
            if mod(iteration,max_iteration)==0
                break;
            end
            if  mod(iteration,5)==0
                if energy_prev<e*1.01
    %                 if constant<5000
    %                     constant=constant*2;
    %                     w=(sum(l1(:).^2+l2(:).^2))/(size(faces,1)*3*2)/constant;
    %                 else
    %                     break
    %                 end
                    break;
                else
                    energy_prev=e;
                end
            end
            if e>2*E
                e
                points_temp=reshape(points_prev,size(points_old1));
                break;
            end
            if alpha==0
                continue
            end
            if e<1e-40
                break
            end

        end

        %     close
        %     figure
        %     plot(time1,log10(e_record));drawnow;
        %     figure
        %     plot(time1,alpha_record);drawnow;
        %     figure
        %     plot(time1,zd_record);drawnow;
        objname=sprintf("SQP\\result\\%d_%d\\%s_%f.obj",energy_type,p,name1(11:end-13),t)
        matname=sprintf("SQP\\result\\%d_%d\\%s_%f.mat",energy_type,p,name1(11:end-13),t);

        save(matname,'points_temp');
        save_obj(points_temp,faces,objname);
        if ceil(t*100)-t*100==0
            figure
            plot(time1,log10(e_record));drawnow;
            figure
            trimesh(faces, points_temp(:,1), points_temp(:,2), points_temp(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');
            drawnow; pause(0.001);
        end
        l_temp=edge_length(faces,points_temp);
        energy=sum((l_temp(:).^2-l_target(:).^2).^2)

    end

end

function grad = compute_grad(points_temp, faces, l_target, energy_type, p,w,initial_points)
    l_temp = edge_length_c(faces, points_temp);
    if energy_type == 1
        grad = compute_grad_c(points_temp, faces, l_target, l_temp);
    elseif energy_type == 2
        grad = sqrt_energy_compute_grad(points_temp, faces, l_target, l_temp);
    elseif energy_type == 3
        grad = energy3_compute_grad(points_temp, faces, l_target);
    elseif energy_type == 4
        grad = energy4_compute_grad(points_temp, faces, l_target);
    elseif energy_type == 5
        weight = 1 ./ l_target.^2;
        grad = energy3_1_compute_grad(points_temp, faces, l_target, weight);      
    elseif energy_type == 6
%         weight = 1 ./ l_target.^p;
%         grad = energyLP_compute_grad(points_temp, faces, l_target, l_temp, weight, p);
        grad = compute_grad_c(points_temp, faces, l_target, l_temp, energy_type, p);
    elseif energy_type==7
        grad = compute_grad_c(points_temp, faces, l_target, l_temp, energy_type, p,w,initial_points);
    end
end

function grad = sqrt_energy_compute_grad(points_temp, faces, l_target, l_temp)
    grad = zeros(size(points_temp));
    for i = 1:size(faces, 1)
        for j = 1:3
            for k = 1:3
                grad(faces(i, j), k) = grad(faces(i, j), k) + sign(l_temp(i, j)^2-l_target(i, j)^2) / sqrt(abs(l_temp(i, j)^2 - l_target(i, j)^2)) * (points_temp(faces(i, j), k) - points_temp(faces(i, mod(j, 3) + 1), k));
                grad(faces(i, mod(j, 3) + 1), k) = grad(faces(i, mod(j, 3) + 1), k) - sign(l_temp(i, j)^2-l_target(i, j)^2) / sqrt(abs(l_temp(i, j)^2 - l_target(i, j)^2)) * (points_temp(faces(i, j), k) - points_temp(faces(i, mod(j, 3) + 1), k));
            end
        end
    end
end

function E = Energy(faces, points, l_target, energy_type, p,w,initial_points)
    l_temp = edge_length_c(faces, points);
    if energy_type == 1
        E = sum((l_temp(:).^2-l_target(:).^2).^2);
    elseif energy_type == 2
        E = sum(abs(l_temp(:).^2 - l_target(:).^2).^(1 / 2));
    elseif energy_type == 3
        E = sum((l_temp(:)-l_target(:)).^2);
    elseif energy_type == 4
        E = sum(exp((l_temp(:)-l_target(:)).^2)-1);
    elseif energy_type == 5
        weight = 1 ./ l_target.^2;
        E = sum(weight(:).*(l_temp(:) - l_target(:)).^2);
    elseif energy_type == 6
        E = sum((l_temp(:)./l_target(:) - 1).^p);
    elseif energy_type==7
        E = sum((l_temp(:).^2-l_target(:).^2).^2);
        E=E+w*sum((initial_points(:)-points(:)).^2);
    end
end

function grad = energy4_compute_grad(points_temp, faces, l_target)
    grad = zeros(size(points_temp));
    l_temp = edge_length(faces, points_temp);
    for i = 1:size(faces, 1)
        for j = 1:3
            for k = 1:3
                g = 2 * exp((l_temp(i, j)-l_target(i, j))^2) * (points_temp(faces(i, j), k) - points_temp(faces(i, mod(j, 3)+1), k)) * (l_temp(i, j) - l_target(i, j)) / l_temp(i, j);
                grad(faces(i, j), k) = grad(faces(i, j), k) + g;
                grad(faces(i, mod(j, 3) + 1), k) = grad(faces(i, mod(j, 3) + 1), k) - g;
            end
        end
    end
end
