function E=energy(faces,point_number,a,points_to_faces,edge_length,point_to_points,A_index)
E=0;
for i=1:point_number
    M=eye(3);
    for j=1:points_to_faces(i,1)      
        theta=0;
        for k=1:3
            if(faces(points_to_faces(i,j+1),k)==i)
                m=points_to_faces(i,j+1);
                theta=acos((edge_length(m,k)^2+edge_length(m,mod(k+1,3)+1)^2-edge_length(m,mod(k,3)+1)^2)/(2*edge_length(m,k)*edge_length(m,mod(k+1,3)+1)));
                break
            end
        end      
        M=M*[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];
        
        theta=a(A_index(i,point_to_points(i,mod(j,points_to_faces(i,1))+2)));
        M=M*[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];
 
    end
    D=M-eye(3);
    E=E+sum(D(:).^2);
    if abs(imag(E))>0
        E
    end
    
end
 