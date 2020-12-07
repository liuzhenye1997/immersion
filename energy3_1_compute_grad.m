function grad=energy3_1_compute_grad(points_temp,faces,l_target,weight)
    grad=zeros(size(points_temp));
    l_temp=edge_length(faces,points_temp);
    for i=1:size(faces,1)
        for j=1:3
            for k=1:3
                grad(faces(i,j),k)=grad(faces(i,j),k)+2*weight(i,j)*(l_temp(i,j)-l_target(i,j))*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k))/l_temp(i,j);
                grad(faces(i,mod(j,3)+1),k)=grad(faces(i,mod(j,3)+1),k)-2*weight(i,j)*(l_temp(i,j)-l_target(i,j))*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k))/l_temp(i,j);
            end
        end
    end
%     grad=compute_grad_c( points_temp, faces, l_target, l_temp);
end