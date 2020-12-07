function points=shape_from_metric(max_iteration,length,faces)
    addpath('splsolver')
    g_is_zero=true;
    points=BunnyFromMetric(max_iteration,length,faces,g_is_zero);
end

