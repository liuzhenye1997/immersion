%local_coordinate为每个点的相邻点在其局部坐标系下的坐标
function local_coordinate=compute_local_coordinate(faces,point_to_points,length,points_to_faces,A)
face_number=size(faces,1);
point_number=size(point_to_points,1);

two_points_to_faces=sparse(faces(:),[faces(:,2);faces(:,3);faces(:,1)],[1:face_number,1:face_number,1:face_number]');

row=zeros(face_number*3,3);
col=zeros(face_number*3,3);
val=zeros(face_number*3,3);

iter=1;
for i=1:point_number
    for j=1:point_to_points(i,1)
        row(iter,:)=[i i i];
        col(iter,:)=3*(point_to_points(i,j+1)-1)+[1 2 3];
        f_id=two_points_to_faces(i,point_to_points(i,j+1));
        for k=1:3
            if(faces(f_id,k)==point_to_points(i,j+1))
                break;
            end
        end
       
        val(iter,1)=length(f_id,mod(k-2,3)+1);
        val(iter,2)=0;
        val(iter,3)=0;       
        iter=iter+1;
    end
    iter=iter-point_to_points(i,1);
    M=eye(3);
    for j=1:point_to_points(i,1)
        val(iter,:)=M*val(iter,:)';
        iter=iter+1;
        theta=0;
        for k=1:3
            if(faces(points_to_faces(i,j+1),k)==i)
                m=points_to_faces(i,j+1);
                theta=acos((length(m,k)^2+length(m,mod(k+1,3)+1)^2-length(m,mod(k,3)+1)^2)/(2*length(m,k)*length(m,mod(k+1,3)+1)));
                break
            end
        end      
        M=M*compute_Z(theta);
        
        theta=A(i,point_to_points(i,mod(j,points_to_faces(i,1))+2));
        M=M*compute_X(theta);
 
    end
end
local_coordinate=sparse(row(:),col(:),val(:));