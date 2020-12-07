%使用中间过程，t从0到1分n份，对迭代终止条件采用wolfe条件，固定一个点。
format long
addpath('splsolver')
% mex  -R2018a sptest5.cpp
max_iteration=1;
name='elephant_S.obj';
% name='temp\\blue_monster_1.000000.obj'
[points_old1,faces,~,~]=readObj(name);
% points_old1=load('blue_monster_S.mat')
% points_old1=points_old1.points_temp;
l1=edge_length(faces,points_old1);
name='elephant_T.obj';
[points_old2,faces,~,~]=readObj(name);
% points_old2(static,:)=points_old1(static,:);
l2=edge_length(faces,points_old2);
energy1_sum=0;
energy2_sum=0;
energy3_sum=0;
logenergy1_sum=0;
logenergy2_sum=0;
logenergy3_sum=0;
error1_max=[];
error2_max=[];
error3_max=[];
energy1=[];
energy2=[];
energy3=[];
start=0;
End=100;
for n=start:End
    t=n*0.01
    l_target=t*l2+(1-t)*l1;
    name1=char(sprintf("output\\max_iteration=100\\1e-10\\2\\elephant_%f.obj",t))
     name2=char(sprintf("output\\max_iteration=100\\1e-10\\2\\elephant_%f.obj",t))
%       name1=char(sprintf("output\\I=1e-0——1e-13\\max_iteration=2000\\u2s_%f.obj",t))
%     name1=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\new_idea\\插值\\LSCM结果\\blue_monster_%f.obj",t));
%     name1=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blue_monster_%f.obj",t));
%     name1=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\LSCM_拟牛顿_阈值1e-10\\blue_monster_%f.obj",t));%1130.2
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\Newton_LSCM_continue极限\\blue_monster_%f.obj",t));%1321.8

%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-5\\blue_monster_%f.obj",t))%1498.4
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-6——1e-7\\blue_monster_%f.obj",t));%1655
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-6——1e-8\\blue_monster_%f.obj",t));%1665.7
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-8\\blue_monster_%f.obj",t))%1555.1
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-4\\blue_monster_%f.obj",t))%1348.9
%     name1=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\max_iteration=2000\\blue_monster_%f.obj",t))%1777.9
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\max_iteration=4000\\blue_monster_%f.obj",t))%1815.7
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\blue_monster_%f.obj",t))%1777.9
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\max_iteration=4000\\blue_monster_%f.obj",t))%1795.8
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-6\\max_iteration=2000\\blue_monster_%f.obj",t))%1637
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\27_3\\max_iteration=2000\\blue_monster_%f.obj",t))
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\27_3\\blue_monster_%f.obj",t))
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-6\\max_iteration=4000\\blue_monster_%f.obj",t))%1690.4
%     name2=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\lu\\blue_monster_%f.obj",t))%17133.2
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_I_Newton_LSCM_continue3\\I=1e-7\\lu_umf\\blue_monster_%f.obj",t))%1763.2
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\LSCM_牛顿_阈值1e-11\\blue_monster_%f.obj",t))
%     name3=char(sprintf("C:\\Users\\liuzhenye\\Desktop\\test2\\拟牛顿法\\2\\blend_liu13_LSCM_拟牛顿_阈值2e-9\\阈值1e-6_max_iteration=2000\\blue_monster_%f.obj",t));
    [points,faces,~,~]=readObj(name1);
    l_temp=edge_length(faces,points);
    energy=sum((l_temp(:).^2-l_target(:).^2).^2)
    energy1_sum=energy1_sum+energy;
    logenergy1_sum=logenergy1_sum+abs(log10(energy));
    energy1=[energy1 energy];
    error1_max=[error1_max max(abs(l_temp(:)./l_target(:)-1))];
    
    [points,faces,~,~]=readObj(name2);
    l_temp=edge_length(faces,points);
    energy=sum((l_temp(:).^2-l_target(:).^2).^2)
    energy2_sum=energy2_sum+energy;
    logenergy2_sum=logenergy2_sum+abs(log10(energy));
    energy2=[energy2 energy];
    error2_max=[error2_max max(abs(l_temp(:)./l_target(:)-1))];
%     
%     [points,faces,~,~]=readObj(name3);
%     l_temp=edge_length(faces,points);
%     energy=sum((l_temp(:).^2-l_target(:).^2).^2)
%     energy3_sum=energy3_sum+energy;
%     logenergy3_sum=logenergy3_sum+abs(log10(energy));
%     energy3=[energy3 energy];
%     error3_max=[error3_max max(abs(l_temp(:)./l_target(:)-1))];
end
plot(1:End-start+1,log10(energy1))
hold on
plot(1:End-start+1,log10(energy2))
hold off
% plot(1:End-start+1,log10(energy3))
% hold off
figure
plot(1:End-start+1,log10(error1_max))
hold on
plot(1:End-start+1,log10(error2_max))
hold off
% plot(1:End-start+1,log10(error3_max))
% hold off
% logerror3_max=sum(abs(log10(error3_max)));
logerror2_max=sum(abs(log10(error2_max)));

logerror1_max=sum(abs(log10(error1_max)));