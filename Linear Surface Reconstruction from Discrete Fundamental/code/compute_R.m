%R则为相邻面坐标系的变换矩阵
function R=compute_R(faces,fr)
face_number=size(faces,1);
R=zeros(face_number,3,3,3);
ad_face=sparse(faces,[faces(size(faces,1)+1:3*size(faces,1)),faces(1:size(faces,1))]',[1:face_number 1:face_number 1:face_number]);
for i=1:face_number
   for j=1:3
        f1=ad_face(faces(i,j),faces(i,mod(j+1,3)+1));
        f2=ad_face(faces(i,mod(j+1,3)+1),faces(i,j));
        R(i,j,:,:)=fr(f1*3-2:3*f1,:)'*fr(f2*3-2:3*f2,:);
   end
end
