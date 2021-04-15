%求能量Ex的系数矩阵，C代表x的系数矩阵,D代表frame的系数矩阵
function [C,D]=conpute_Ex_matrix(point_number,faces,a)
face_number=size(faces,1);

a1=a(:,1:3);
a2=a(:,4:6);

col=zeros(face_number,3,2);
row=zeros(face_number,3,2);
val=zeros(face_number,3,2);
for i=1:face_number
    for j=1:3
        col(i,j,:)=3*(i-1)+j;
        row(i,j,1)=faces(i,mod(j+1,3)+1);
        row(i,j,2)=faces(i,j);
        val(i,j,1)=1;
        val(i,j,2)=-1;
    end
end
C=sparse(row(:),col(:),val(:),point_number,3*face_number);

col=zeros(face_number,3,2);
row=zeros(face_number,3,2);
val=zeros(face_number,3,2);
for i=1:face_number
    for j=1:3
        col(i,j,:)=3*(i-1)+j;
        row(i,j,1)=3*(i-1)+1;
        row(i,j,2)=3*(i-1)+2;
        val(i,j,1)=a1(i,j);
        val(i,j,2)=a2(i,j);
    end
end
D=sparse(row(:),col(:),val(:),3*face_number,3*face_number);