%Mk为论文中的Mk
function Mk=compute_Mk(faces,point_number,A,points_to_faces,edge_length,point_to_points)
Mk=zeros(point_number*3,3*size(point_to_points,2));
face_number=size(faces,1);

for i=1:point_number
    M=eye(3);
    X=zeros(points_to_faces(i,1)*3,3);
    X(3*(1:points_to_faces(i,1))-2,1)=1;
    j=1:points_to_faces(i,1);
    theta=A(i,point_to_points(i,mod(j,points_to_faces(i,1))+2));
    X(3*j-1,2)=cos(theta);
    X(3*j,3)=cos(theta);
    X(3*j-1,3)=-sin(theta);
    X(3*j,2)=sin(theta);
    Z=zeros(points_to_faces(i,1)*3,3);
    for j=1:points_to_faces(i,1)      
        theta=0;
        for k=1:3
            if(faces(points_to_faces(i,j+1),k)==i)
                m=points_to_faces(i,j+1);
                theta=acos((edge_length(m,k)^2+edge_length(m,mod(k+1,3)+1)^2-edge_length(m,mod(k,3)+1)^2)/(2*edge_length(m,k)*edge_length(m,mod(k+1,3)+1)));
                break
            end
        end      
        Z(3*j-2:3*j,:)=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
        M=M*Z(3*j-2:3*j,:);
        if j~=points_to_faces(i,1)      
            M=M*X(3*j-2:3*j,:);           
        end
        
    end
    Mk(3*i-2:3*i,1:3)=M;
    for j=1:points_to_faces(i,1)-1
        index=mod(j+points_to_faces(i,1)-2,points_to_faces(i,1))+1;
        M=M*X(3*index-2:3*index,:); 
        M=M*Z(3*j-2:3*j,:);
        M=Z(3*j-2:3*j,:)'*M;        
        M=X(3*j-2:3*j,:)'*M; 
        
        index=j+1;
        Mk(3*i-2:3*i,3*(index)-2:3*(index))=M;
    end

end

%将Mk转换成稀疏矩阵形式，方便检索
col=zeros(face_number*3,3,3);
row=zeros(face_number*3,3,3);
val=zeros(face_number*3,3,3);
iter=1;
for i=1:point_number
    for j=1:points_to_faces(i,1)    
        row(iter,:,:)=3*(i-1)+[1 1 1;2 2 2;3 3 3];
        col(iter,:,:)=3*(point_to_points(i,j+1)-1)+[1 2 3;1 2 3;1 2 3];
        val(iter,:,:)=Mk(3*i-2:3*i,3*j-2:3*j);
        iter=iter+1;
    end
end
Mk=sparse(row(:),col(:),val(:));


