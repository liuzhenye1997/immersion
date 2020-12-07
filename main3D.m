%目前只能应用在三维封闭网格上
addpath('shapefrommetric')
max_iteration=80;
target_name=char("bunny.obj");
% % Output_path="output\\I=1e-0——1e-13\\spirals";
Output_path="output\\bunny";
show_energy=true;%是否展示迭代过程中的能量变化图
show_result=true;%s是否展示结果图
energy_type=5;%第一种能量是四次能量，第二种能量是1/4次能量，第三种能量是二次能量，即最近一次的
initialization_obj="BunnyFromMetric_60.obj";%初始化使用的网格，可以不使用
immersion3D(max_iteration,target_name,Output_path,show_energy,show_result,energy_type);
