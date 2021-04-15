%a1,a2,a3��Ϊ�����ε���������fr1�������ϵ�ϵ�ϵ����a3��Ϊ0
function [a1,a2]=compute_a(points,faces,f)
a1=zeros(size(faces));
a2=zeros(size(faces));
a3=zeros(size(faces));
for i=1:size(faces,1)
    for j=1:3
        edge=points(faces(i,mod(j+1,3)+1),:)-points(faces(i,j),:);  
        a1(i,j)=dot(f(3*i-2:3*i,1),edge);
        a2(i,j)=dot(f(3*i-2:3*i,2),edge);
        a3(i,j)=dot(f(3*i-2:3*i,3),edge);
    end
end

