%�Զ�ά������в�ֵ
max_iteration = 100;
source_name = char("spirals_S.obj");
target_name = char("spirals_T.obj");
% % Output_path="output\\I=1e-0����1e-13\\spirals";
Output_path = "output\\spirals";
frames = 100;
show_energy = true; %�Ƿ�չʾ���������е������仯ͼ
show_result = true; %s�Ƿ�չʾ���ͼ
energy_type = 6; %��һ���������Ĵ��������ڶ���������1/4�������������������Ƕ�������,���ģ��������������Կ���L2������������������LP����
p = 2;

immersion(max_iteration, source_name, target_name, frames, Output_path, show_energy, show_result, energy_type, p);