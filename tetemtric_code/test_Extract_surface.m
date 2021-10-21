[x1, tets] = readMESH( basedir + "source.mesh" );
[p,t]=Extract_surface(x1,tets);
% for i=1:size(p,1)
%     t=[t ;[i i i]];
% end
T0=triangulation(t,p(:,1),p(:,2),p(:,3));
trimesh(T0)
drawnow