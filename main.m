max_iteration=100;
source_name=char("spirals_S.obj");
target_name=char("spirals_T.obj");
% % Output_path="output\\I=1e-0����1e-13\\spirals";
Output_path="output\\spirals";
frames=100;
show_energy=false;%�Ƿ�չʾ���������е������仯ͼ
show_result=false;%s�Ƿ�չʾ���ͼ

immersion(max_iteration,source_name,target_name,frames,Output_path,show_energy,show_result);