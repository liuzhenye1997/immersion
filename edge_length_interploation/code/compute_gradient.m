%gradient为能量的梯度，gradient(i,j)则为A(i,j)对应的梯度
function gradient=compute_gradient(Mk,faces,A)
face_number=size(faces,1);
row=zeros(face_number,3,2);
col=zeros(face_number,3,2);
val=zeros(face_number,3,2);
point_number=max(faces(:));
for i=1:face_number
    for j=1:3
        row(i,j,1)=faces(i,j);
        row(i,j,2)=faces(i,mod(j,3)+1);
        col(i,j,2)=faces(i,j);
        col(i,j,1)=faces(i,mod(j,3)+1);
%         theta=A(faces(i,j),faces(i,mod(j,3)+1));
%         val(i,j,:)=(Mk(3*faces(i,j),3*faces(i,mod(j,3)+1)-1)-Mk(3*faces(i,j)-1,3*faces(i,mod(j,3)+1)))*cos(theta)+(Mk(3*faces(i,j),3*faces(i,mod(j,3)+1))+Mk(3*faces(i,j)-1,3*faces(i,mod(j,3)+1)-1))*sin(theta);
    end
end
i=1:face_number;
j=1:3;
r1=3*faces(i,j);
r2=3*faces(i,j)-1;
c1=3*faces(i,mod(j,3)+1)-1;
c2=3*faces(i,mod(j,3)+1);
index1=sub2ind([point_number*3 point_number*3],r1,c1);
index2=sub2ind([point_number*3 point_number*3],r2,c2);
index3=sub2ind([point_number*3 point_number*3],r1,c2);
index4=sub2ind([point_number*3 point_number*3],r2,c1);
theta=A(sub2ind([point_number point_number],faces(i,j),faces(i,mod(j,3)+1)));
val1=zeros(face_number,3,2);
val1(:,:,1)=(Mk(index1)-Mk(index2)).*cos(theta)+(Mk(index3)+Mk(index4)).*sin(theta);
val1(:,:,2)=val1(:,:,1);
gradient=sparse(row(:),col(:),val1(:));