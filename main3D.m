%Ŀǰֻ��Ӧ������ά���������
addpath('shapefrommetric')
max_iteration=80;
target_name=char("bunny.obj");
% % Output_path="output\\I=1e-0����1e-13\\spirals";
Output_path="output\\bunny";
show_energy=true;%�Ƿ�չʾ���������е������仯ͼ
show_result=true;%s�Ƿ�չʾ���ͼ
energy_type=5;%��һ���������Ĵ��������ڶ���������1/4�������������������Ƕ��������������һ�ε�
initialization_obj="BunnyFromMetric_60.obj";%��ʼ��ʹ�õ����񣬿��Բ�ʹ��
immersion3D(max_iteration,target_name,Output_path,show_energy,show_result,energy_type);