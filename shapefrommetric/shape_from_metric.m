%����shape_from_metric����������߳���һ������Ƕ��
function points = shape_from_metric(max_iteration, length, faces)
    addpath('splsolver')
    %����ԭ���ߴ������������Ĳ�����g�ڴ�����һֱΪ0������ֱ����Ϊ0���Խ�ʡʱ��
    g_is_zero = true;
    points = BunnyFromMetric(max_iteration, length, faces, g_is_zero);
end
