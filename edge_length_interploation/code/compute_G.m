%%G为每个点的坐标架在世界坐标系下的值
function G=compute_G(G1,M,faces)
face_number=size(faces,1);
row=zeros(face_number,3,3,3);
col=zeros(face_number,3,3,3);
val=zeros(face_number,3,3,3);
for i=1:face_number
    for j=1:3
        row(i,j,:,:)=3*(3*(i-1)+j-1)+[1 1 1;2 2 2;3 3 3];
        col(i,j,:,:)=3*faces(i,mod(j,3)+1)-3+[1 2 3;1 2 3;1 2 3];
        val(i,j,:,:)=M(3*faces(i,j)-2:3*faces(i,j),3*faces(i,mod(j,3)+1)-2:3*faces(i,mod(j,3)+1));
    end
end

A=sparse(row(:),col(:),val(:));

for i=1:face_number
    for j=1:3
        row(i,j,:,:)=3*(3*(i-1)+j-1)+[1 1 1;2 2 2;3 3 3];
        col(i,j,:,:)=3*faces(i,j)-3+[1 2 3;1 2 3;1 2 3];
        val(i,j,:,:)=eye(3);
    end
end
B=sparse(row(:),col(:),val(:));
A=A-B;

A1=A(:,1:3);
A2=A(:,4:end);
G2=-(A2'*A2)\(A2'*A1*G1);
G=[G1;G2];
d=A*G;
sum(abs(d(:)));

for i=2:(size(G,1)/3)
    G(3*i-2:3*i,:)=G(3*i-2:3*i,:)';
end





