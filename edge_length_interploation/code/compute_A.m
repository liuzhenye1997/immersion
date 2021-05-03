%A为二面角，A(i,j)即为点i和点j的连线对应的二面角
function A=compute_A(faces,points)
face_number=size(faces,1);

%two_points_to_points指的是找到三角形的两点后得到第三个点
two_points_to_points=sparse(faces(:),[faces(:,2);faces(:,3);faces(:,1)],[faces(:,3);faces(:,1);faces(:,2)]);

col=zeros(face_number,3);
row=zeros(face_number,3);
val=zeros(face_number,3);
for i=1:face_number
    for j=1:3
        row(i,j)=faces(i,j);
        col(i,j)=faces(i,mod(j,3)+1);
        a=points(two_points_to_points(faces(i,mod(j,3)+1),faces(i,j)),:);
        b=points(faces(i,mod(j,3)+1),:);
        c=points(faces(i,j),:);
        d=points(two_points_to_points(faces(i,j),faces(i,mod(j,3)+1)),:);
        val(i,j)=dihedral_angle(a,b,c,d);
    end
end
val(abs(imag(val))<1e-5)=real(val);
A=sparse(row(:),col(:),(val(:)));