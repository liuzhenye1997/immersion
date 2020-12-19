%%利用shape_from_metric初始化,然后再选择某一种能量进行优化
function immersion3D(max_iteration, target_name, Output_path, show_energy, show_result, energy_type, p, initialization_obj)
    format long
    addpath('splsolver')
    if energy_type ~= 1 && energy_type ~= 2 && energy_type ~= 3 && energy_type ~= 4 && energy_type ~= 5 && energy_type ~= 6
        error("energy_type应为1,2,3,4,5或者6");
    end
    if nargin == 6
        p = 2;
    end
    [points_old1, faces, ~, ~] = readObj(target_name);
    
    l1 = edge_length(faces, points_old1);
    
    point_number = size(points_old1, 1);
    
    fix = [1, point_number + 1, point_number * 2 + 1];
    free = 1:3 * point_number;
    free(fix) = [];
    
    l_target = l1;
    time = [];
    E_record = [];
    t_record = [];
    Time = [];
    Alpha_record = [];
    ZD_record = [];
    
    T1 = clock;
    threshold = 1e-19;
  
    for t = 0
        l_target = l1;
        if nargin == 7 || nargin == 6
            points_temp = shape_from_metric(10, l_target, faces);
            
        elseif nargin == 8
            [points_temp, ~] = readObj(char(initialization_obj));
        end
        
        H = eye(size(points_old1, 1)*3);
        I = sparse(1:size(points_old1, 1)*3, 1:size(points_old1, 1)*3, 1e-6);
        energy_prev = Energy(faces, points_temp, l_target, energy_type, p)
        

        e_record = [];
        e_record = [e_record, energy_prev];
        alpha_record = [];
        alpha_record = [alpha_record, 0];
        time = [time, 0];
        t1 = clock;
        time1 = 0;
        iteration = 0;
        
        while 1
            iteration = iteration + 1;
            if iteration == 1
                if energy_type == 5 || energy_type == 4 || energy_type == 3 || energy_type == 2 || energy_type == 1
                    [row, col, ~, H] = svd_positive_hessian_matrix_c(points_temp, faces, l_target, energy_type);
                elseif energy_type == 6
                    [row, col, ~, H] = svd_positive_hessian_matrix_c(points_temp, faces, l_target, energy_type, p);
                end
                H = H + sparse(1:size(points_old1, 1)*3, 1:size(points_old1, 1)*3, 10^-15);
                H_index = ij2nzIdxs(H, uint64(row(:)), uint64(col(:)));
                [M, N] = find(H);
                k = unique([find(sum(M(:) == fix, 2)), find(sum(N(:) == fix, 2))]);
                H1_index = 1:size(M, 1);
                H1_index(k) = [];
                H1 = H(free, free);
            else
                if energy_type == 5 || energy_type == 4 || energy_type == 3 || energy_type == 2 || energy_type == 1
                    value = svd_positive_hessian_matrix_value_c(points_temp, faces, l_target, energy_type);
                elseif energy_type == 6
                    value = svd_positive_hessian_matrix_value_c(points_temp, faces, l_target, energy_type, p);
                end
                H_val = myaccumarray(H_index, value);
                update_matrix_single_thread(H1, H_val, H1_index)
                
            end
            if iteration == 1
                if energy_type ~= 2
                    solver = splsolver(H1, 'ldlt');
                else
                    solver = splsolver(H1, 'ldlt');
                end
            end
            if energy_type ~= 6
                grad = compute_grad(points_temp, faces, l_target, energy_type);
            else
                grad = compute_grad(points_temp, faces, l_target, energy_type, p);
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
            E = Energy(faces, points_temp, l_target, energy_type, p)
            if energy_type == 1
                alpha = best_step(points_temp, faces, d, l_target);
            else
                alpha = 100;
                for i = 1:100
                    points_temp1 = points_temp;
                    points_temp1(2:end, :) = points_temp1(2:end, :) - alpha * d(2:end, :);
                    if Energy(faces, points_temp1, l_target, energy_type, p) < E
                        break;
                    else
                        alpha = alpha / 2;
                    end
                end
            end
            
            points_prev = points_temp;
            points_temp(free) = points_temp(free) - alpha * d(free);
            points_temp = reshape(points_temp, size(points_old1));
            
            
            e = Energy(faces, points_temp, l_target, energy_type, p);
            e_record = [e_record, e];
            time1 = [time1, etime(clock, t1)];
            alpha_record = [alpha_record, alpha];
            
            if mod(iteration, max_iteration) == 0
                break;
            end
            if mod(iteration, 5) == 0
                if energy_prev < e * 1.01
                    break;
                else
                    energy_prev = e;
                end
            end
            if e > 2 * E
                e
                points_temp = reshape(points_prev, size(points_old1));
                break;
            end
            if alpha == 0
                break;
            end
            if e < 1e-40
                break
            end
            
        end
 
        if show_energy
            figure
            plot(time1, log10(e_record));
            drawnow;
        end
        objname = Output_path + sprintf("_%d.obj", max_iteration);
        matname = Output_path + sprintf("_%d.mat", max_iteration);
        save(matname, 'points_temp');
        save_obj(points_temp, faces, objname);
        if show_result
            figure
            trimesh(faces, points_temp(:, 1), points_temp(:, 2), points_temp(:, 3), 'edgecolor', 'k');
            axis off;
            axis equal;
            title('output');
            drawnow;
            pause(0.001);
        end
    end
    
    mesh = Triangular_mesh(faces, l1);
    [max_abs_rel_err, min_rel_err, max_rel_err, l2_rel_err, max_log_err, min_log_err] = evaluate(mesh.f_number, mesh.v_dst, mesh.v_src, points_temp, mesh.length);
    fprintf("最大相对误差为%f,L2误差为%f", max_abs_rel_err, l2_rel_err);
end

%%
%计算各种能量对应的梯度
function grad = compute_grad(points_temp, faces, l_target, energy_type, p)
    l_temp = edge_length_c(faces, points_temp);
    if energy_type == 5 || energy_type == 4 || energy_type == 3 || energy_type == 2 || energy_type == 1
        grad = compute_grad_c(points_temp, faces, l_target, l_temp, energy_type);
    elseif energy_type == 6
        grad = compute_grad_c(points_temp, faces, l_target, l_temp, energy_type, p);
    end
end

%%
%计算各种能量
function E = Energy(faces, points, l_target, energy_type, p)
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
        E = sum((l_temp(:)./l_target(:)-1).^p);
    end
end
