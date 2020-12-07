 name='bunny.obj';
[points,faces,~,~]=readObj(name);
l_target=edge_length(faces,points); 
name=sprintf("bunny_5000.obj",t);
[points,faces,~,~]=readObj(char(name));

grad=compute_grad(points,faces,l_target);
sprintf("梯度绝对值的平均值为%e\n",sum(abs(grad(:)))/(3*size(grad,1)))


function grad=compute_grad(points_temp,faces,l_target)
    grad=zeros(size(points_temp));
    l_temp=edge_length(faces,points_temp);
    for i=1:size(faces,1)
        for j=1:3
            for k=1:3
                grad(faces(i,j),k)=grad(faces(i,j),k)+2*(l_temp(i,j)^2-l_target(i,j)^2)*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
                grad(faces(i,mod(j,3)+1),k)=grad(faces(i,mod(j,3)+1),k)-2*(l_temp(i,j)^2-l_target(i,j)^2)*(points_temp(faces(i,j),k)-points_temp(faces(i,mod(j,3)+1),k));
            end
        end
    end
end