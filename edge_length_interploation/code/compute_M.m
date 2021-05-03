%M为相邻两点的局部坐标系的变换矩阵
function M=compute_M(A,faces,point_to_points,edge_length,points_to_faces)
face_number=size(faces,1);
point_number=max(faces(:));
min_id=ones(point_number,1).*(point_number+1);
for i=1:face_number
    for j=1:3
        min_id(faces(i,j))=min(min_id(faces(i,j)),faces(i,mod(j,3)+1));
        min_id(faces(i,j))=min(min_id(faces(i,j)),faces(i,mod(j-2,3)+1));
    end
end

row=zeros(face_number,3,3,3);
col=zeros(face_number,3,3,3);
val=zeros(face_number,3,3,3);
for i=1:face_number
    for j=1:3
        src=faces(i,j);
        dst=faces(i,mod(j,3)+1);
        row(i,j,:,:)=3*(src-1)+[1 1 1;2 2 2;3 3 3];
        col(i,j,:,:)=3*(dst-1)+[1 2 3;1 2 3;1 2 3];
        n=2;
        M=eye(3);
        while 1
            
            theta=0;
            for k=1:3
                if(faces(points_to_faces(src,n),k)==faces(i,j))
                    m=points_to_faces(src,n);
                    theta=acos((edge_length(m,k)^2+edge_length(m,mod(k+1,3)+1)^2-edge_length(m,mod(k,3)+1)^2)/(2*edge_length(m,k)*edge_length(m,mod(k+1,3)+1)));
                    break
                end
            end      
            M=M*[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
            n=n+1;
            if (point_to_points(src,mod(n-2,point_to_points(src,1))+2)~=dst)
                theta=A(src,point_to_points(src,mod(n-2,point_to_points(src,1))+2));
                M=M*[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
            else
                break;
            end
            
        end
        M=M*compute_Z(pi);
        index=2;
        while(point_to_points(faces(i,mod(j,3)+1),index)~=faces(i,j))
            index=index+1;
            index=mod(index-2,points_to_faces(faces(i,mod(j,3)+1),1))+2;
        end
        if index>2
            for n=index:point_to_points(faces(i,mod(j,3)+1),1)+1
                theta=0;
                for k=1:3
                    if(faces(points_to_faces(faces(i,mod(j,3)+1),n),k)==faces(i,mod(j,3)+1))
                        m=points_to_faces(faces(i,mod(j,3)+1),n);
                        theta=acos((edge_length(m,k)^2+edge_length(m,mod(k+1,3)+1)^2-edge_length(m,mod(k,3)+1)^2)/(2*edge_length(m,k)*edge_length(m,mod(k+1,3)+1)));
                        break
                    end
                end    
                M=M*[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];

                theta=A(faces(i,mod(j,3)+1),point_to_points(faces(i,mod(j,3)+1),mod(n-1,points_to_faces(faces(i,mod(j,3)+1),1))+2));
                M=M*[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
            end
        end
        val(i,j,:,:)=M;
    end
end
M=sparse(row(:),col(:),val(:));

