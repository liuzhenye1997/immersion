%利用shape_from_metric方法求给定边长的一个近似嵌入
function points = shape_from_metric(max_iteration, length, faces)
    addpath('splsolver')
    %由于原作者代码有误，与论文不符，g在代码中一直为0，所以直接设为0，以节省时间
    g_is_zero = true;
    points = BunnyFromMetric(max_iteration, length, faces, g_is_zero);
end
