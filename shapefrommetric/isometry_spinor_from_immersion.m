%论文5.4.1,算法8.此处设置了a=0。
%将仿射矩阵极分解,提取出旋转矩阵,然后求出对应的lambda
function lambda = isometry_spinor_from_immersion(dt, length, face_number, vertex_theta, vertex_src, vertex_dst, vertex_next, face_hedge, points, faces, lambda)
    %求仿射矩阵df,利用原网格和目标网格的法向和两条边求
    normal = compute_normal(points, faces);
    hedge0 = face_hedge';
    vertex_next_temp = vertex_next';
    hedge1 = vertex_next_temp(hedge0);
    vertex_src_temp = vertex_src';
    vertex_dst_temp = vertex_dst';
    P_hedge0src = points(vertex_src_temp(hedge0), :);
    P_hedge0dst = points(vertex_dst_temp(hedge0), :);
    P_hedge1src = points(vertex_src_temp(hedge1), :);
    P_hedge1dst = points(vertex_dst_temp(hedge1), :);
    V0 = P_hedge0dst - P_hedge0src;
    V1 = P_hedge1dst - P_hedge1src;
    

    l0 = length(hedge0);
    l1 = length(hedge1);
    length_V0 = sqrt(sum(V0.^2, 2));
    
    %N为法向,E为一条边的单位方向，NE为那条边上的垂线
    N = normal;
    E = V0 ./ length_V0;
    NE = [N(:, 2) .* E(:, 3) - N(:, 3) .* E(:, 2), N(:, 3) .* E(:, 1) - N(:, 1) .* E(:, 3), N(:, 1) .* E(:, 2) - N(:, 2) .* E(:, 1)];
    u = sum(E.*V1, 2);
    v = sum(NE.*V1, 2);
    vertex_theta_temp = vertex_theta';
    theta = vertex_theta_temp(hedge1);
    
    a = (u - l1 ./ l0 .* length_V0 .* cos(theta)) ./ (l1 .* sin(theta));
    b = v ./ (l1 .* sin(theta));
    F = a .* E + b .* NE;
    
    %df为每个面对应的变换矩阵。即该矩阵能将原本在平面上的三角形面片映射到现在所在的位置上
    df = zeros(face_number*3, 3);
    df(1:3:face_number*3, :) = N;
    df((1:3:face_number*3)+1, :) = V0 ./ l0;
    df((1:3:face_number*3)+2, :) = F;
    
    %对矩阵df进行极分解,U为分解后得到的正交矩阵。然后再变换成四元数
    U = polar_decomposition(df, face_number, F, V0, length_V0, l0);
    q = Quaternion.batch_matrix_to_q(U);
    
    %由于开方存在正负两种情况，故得到的lambda也存在正负两种情况，故让其与刚开始的lambda的内积大于0
    prod = lambda(:, 1) .* q(:, 1) + lambda(:, 2) .* q(:, 2) + lambda(:, 3) .* q(:, 3) + lambda(:, 4) .* q(:, 4);
    for i = 1:face_number
        if prod(i) < 0
            q(i, :) = -q(i, :);
        end
    end
    lambda_pre = lambda;
    lambda = q;
    
    %由于a=0，故事实上下面并没有改变lambda;
    a = 0;
    lambda = penalty_parameter(lambda_pre, lambda, a, dt);
end

%% 算法8第7,8行
function lambda = penalty_parameter(lambda_pre, lambda, a, dt)
    if abs(a) < 1 * exp(-6)
        r = 0;
    else
        r = exp(-dt/a);
    end
    
    lambda = lambda * (1 - r);
    lambda = lambda + r * lambda_pre;
    lambda_norm = sqrt(sum(lambda.^2, 2));
    lambda = lambda ./ lambda_norm;
end

%% 对矩阵df进行极分解,U为分解后得到的正交矩阵。
%由于是对该矩阵专门进行的优化，利用了它的结构的特殊性，故看起来比较复杂。
%假设df=UP,U是正交阵，P为半正定矩阵。则P=sqrt(df'*df),U=A*P_inv,P_inv为P的逆矩阵。
%因此该函数先求ddf=df'*df，然后开方得到P，再求得P_inv，进而得到U。df的极分解完成
%由于单位法向量映射到单位法向量，故ddf=df'*df形如
%[1 0 0
% 0 a b
% 0 b c].
%因此实际上是对2*2的矩阵进行开方得到P_2
function U = polar_decomposition(df, face_number, F, V0, length_V0, l0)   
    ddf = zeros(face_number, 4);
    ddf(:, 1) = (length_V0 ./ l0).^2;
    ddf(:, 4) = sum(F.^2, 2);
    ddf(:, 3) = sum(F.*V0, 2) ./ l0;
    ddf(:, 2) = (sum(F .* V0, 2)) ./ l0;
    
    P_inv = zeros(3*face_number, 3);
    P_inv(1:3:3*face_number, 1) = 1;
    P_2 = batch_sqrt2_2(ddf);
    
    P_det = P_2(:, 1) .* P_2(:, 4) - P_2(:, 3).^2;
    P_inv((1:3:3*face_number)+1, 2) = P_2(:, 4) ./ P_det;
    P_inv((1:3:3*face_number)+1, 3) = -P_2(:, 3) ./ P_det;
    P_inv((1:3:3*face_number)+2, 2) = -P_2(:, 2) ./ P_det;
    P_inv((1:3:3*face_number)+2, 3) = P_2(:, 1) ./ P_det;
    
    U = zeros(face_number, 9);
    A = zeros(face_number, 9);
    B = reshape(P_inv', 9, face_number)';
    for i = 1:3
        for j = 1:3
            A(:, 3*i+j-3) = df((0:3:3*face_number-1)+j, i);
        end
    end
    for i = 1:3
        for j = 1:3
            U(:, 3*i+j-3) = sum(A(:, 3 * i - 2:3 * i).*B(:, 3 * j - 2:3 * j), 2);
        end
    end
end

%% 对2*2的对角阵同时进行开方
function result = batch_sqrt2_2(A)
    sin_theta = A(:, 3) ./ sqrt(A(:, 1).*A(:, 4));
    cos_theta = sqrt(1-sin_theta.^2);
    c = sqrt(A(:, 1)+A(:, 4)+2*sqrt(A(:, 1) .* A(:, 4)).*cos_theta);
    h = A(:, 3) ./ c;
    a = sqrt(A(:, 1)-h.^2);
    b = sqrt(A(:, 4)-h.^2);
    result = [a, h, h, b];
end
