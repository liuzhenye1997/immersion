%A为矩阵形式，A(i,j)为点i和点j连线对应的二面角，方便检索
%a为向量形式，即将二面角排成一列，方便运算
%a(A_index)=A;
function [a,A_index]=A_to_a(faces,A)
face_number=size(faces,1);
a=zeros(3*face_number/2,1);
iter=1;
row=zeros(3*face_number/2,2);
col=zeros(3*face_number/2,2);
val=zeros(3*face_number/2,2);
for i=1:face_number
    for j=1:3
        if A(faces(i,j),faces(i,mod(j,3)+1))~=100
            a(iter)=A(faces(i,j),faces(i,mod(j,3)+1));
            A(faces(i,j),faces(i,mod(j,3)+1))=100;
            A(faces(i,mod(j,3)+1),faces(i,j))=100;
            row(iter,1)=faces(i,j);
            col(iter,1)=faces(i,mod(j,3)+1);
            row(iter,2)=faces(i,mod(j,3)+1);
            col(iter,2)=faces(i,j);
            val(iter,:)=iter;
            
            iter=iter+1;
        end
    end
end
A_index=sparse(row(:),col(:),val(:));
