%重建网格
function points=reconstruction(faces,local_coordinate,G,x0)
face_number=size(faces,1);
row=zeros(face_number,3,2);
col=zeros(face_number,3,2);
val=zeros(face_number,3,2);
for i=1:face_number
    for j=1:3
        row(i,j,:)=3*(i-1)+j;
        col(i,j,1)=faces(i,j);
        col(i,j,2)=faces(i,mod(j,3)+1);
        val(i,j,1)=1;
        val(i,j,2)=-1;
    end
end
A=sparse(row(:),col(:),val(:));

b=zeros(3*face_number,3);
for i=1:face_number
    for j=1:3
%         for k=1:point_to_points(faces(i,mod(j,3)+1),1)
%             if(point_to_points(faces(i,mod(j,3)+1),k+1)==faces(i,j))
%                 break;
%             end
%         end
        b(3*(i-1)+j,:)=G(3*faces(i,mod(j,3)+1)-2:3*faces(i,mod(j,3)+1),:)*local_coordinate(faces(i,mod(j,3)+1),3*faces(i,j)-2:3*faces(i,j))';
    end
end

% point_number=size(A,1);
% A1=A(:,index_static);
% A2=A(:,index_free);
% x1=(A2'*A2)\(A2'*(b-A1*x0));
% points=zeros(point_number,3);
% points(index_free,:)=x1;
% points(index_static,:)=x0;

A1=A(:,1);
A2=A(:,2:end);
x1=(A2'*A2)\(A2'*(b-A1*x0));
points=[x0;x1];
