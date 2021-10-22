format long
addpath('splsolver')

energy_type=7;
p=1;
list=dir("./tetinterp_meshes");
path=sprintf('./SQP/result/%d_%d',energy_type,p);
if ~exist(path,'dir')
    mkdir(path);
end
path=sprintf('./SQP/obj/%d_%d',energy_type,p);
if ~exist(path,'dir')
    mkdir(path);
end
d=size(list,1)/2;
for obj=6:6
    [~,faces]=tetmetric(list(2*obj-1).name(1:end-12));
    name1=sprintf('./SQP/obj/%s_0.000000.obj',list(2*obj-1).name(1:end-12));

    name2=sprintf('./SQP/obj/%s_1.000000.obj',list(2*obj-1).name(1:end-12));
    interpolation(name1,name2,energy_type,p);
end




