%计算原网格的每个面上的坐标系框架，即fr
function fr=compute_frame(points,faces)
faces_number=size(faces,1);
fr=zeros(3*faces_number,3);
length=edge_length(faces,points);
%以第一条边为x轴，改编的垂线为y轴，法线为z轴
for i=1:faces_number
    b=points(faces(i,3),:)-points(faces(i,1),:);
    a=points(faces(i,2),:)-points(faces(i,1),:);
    fr(3*(i-1)+1:3*i,1)=a./length(i,1);
    dot=sum(b.*a);
    dot=sum(dot(:));
    h=b-dot*a./(length(i,1).^2);
    h=h./norm(h);
    fr(3*(i-1)+1:3*i,2)=h;
    fr(3*(i-1)+1:3*i,3)=cross(fr(3*(i-1)+1:3*i,1),fr(3*(i-1)+1:3*i,2));
end
