%对二维网格进行插值
max_iteration = 100;
source_name = char("spirals_S.obj");
target_name = char("spirals_T.obj");
% % Output_path="output\\I=1e-0——1e-13\\spirals";
Output_path = "output\\spirals";
frames = 100;
show_energy = true; %是否展示迭代过程中的能量变化图
show_result = true; %s是否展示结果图
energy_type = 6; %第一种能量是四次能量，第二种能量是1/4次能量，第三种能量是二次能量,第四，五种能量都可以看做L2能量，第六种能量是LP能量
p = 2;

immersion(max_iteration, source_name, target_name, frames, Output_path, show_energy, show_result, energy_type, p);