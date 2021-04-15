function show_frame(points,faces,fr)
mid=zeros(size(faces,1),3);
for i=1:size(faces,1)
    mid(i,:)=sum(points(faces(i,:),:))/3;
    
end
points_old1=points;
mid=zeros(size(faces,1),3);
for i=1:size(faces,1)
    mid(i,:)=sum(points_old1(faces(i,:),:))/3;   
end
fr=fr*0.2;
for i=1:3
    plot3([mid(:,1),mid(:,1)+fr(1:3:end,i)]',[mid(:,2),mid(:,2)+fr(2:3:end,i)]',[mid(:,3),mid(:,3)+fr(3:3:end,i)]');
    hold on
end
% plot3([mid(:,1),points_old1(faces(:,1),1)]',[mid(:,2),points_old1(faces(:,1),2)]',[mid(:,3),points_old1(faces(:,1),3)]');
% hold on
% plot3([mid(:,1),points_old1(faces(:,1),1)]',[mid(:,2),points_old1(faces(:,1),2)]',[mid(:,3),points_old1(faces(:,1),3)]');
% hold on
trimesh(faces,points(:,1),points(:,2),points(:,3))
axis equal
axis off
hold off