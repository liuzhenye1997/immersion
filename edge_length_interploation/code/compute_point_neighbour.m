%计算每个点相邻的面和每个点相邻的面
function [point_to_faces,point_to_points]=compute_point_neighbour(faces)
degree = points_degree(faces);
point_number = max(faces(:));
max_degree=max(degree(:));
face_number=size(faces,1);

%two_points_to_points指的是找到三角形的两点后得到第三个点
%two_points_to_faces指的是找到三角形的两点后得到这个面的序号
two_points_to_points=sparse(faces(:),[faces(:,2);faces(:,3);faces(:,1)],[faces(:,3);faces(:,1);faces(:,2)]);
two_points_to_faces=sparse(faces(:),[faces(:,2);faces(:,3);faces(:,1)],[1:face_number,1:face_number,1:face_number]');
point_to_faces=ones(point_number,max_degree+1).*(point_number+1);
point_to_faces(:,1)=degree;
%寻找点i的一个邻点作为初始点
for i=1:face_number
    for j=1:3
        point_to_faces(faces(i,j),2)=min(point_to_faces(faces(i,j),2),faces(i,mod(j,3)+1));
    end
end

for i=1:point_number
    for j=1:degree(i)-1
        point_to_faces(i,2+j)=two_points_to_points(i,point_to_faces(i,1+j));
    end
end
point_to_points=point_to_faces;
for i=1:point_number
    for j=1:degree(i)
        point_to_faces(i,1+j)=two_points_to_faces(i,point_to_faces(i,1+j));
    end
end


