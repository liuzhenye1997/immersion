function [p,t]=Extract_surface(x,tets)
pn=size(x,1);
f1=tets(:,1:3);
f2=tets(:,2:4);
f3=[tets(:,3:4) tets(:,1)];
f4=[tets(:,4) tets(:,1:2)];
faces=[f1;f2;f3;f4];
faces_temp=sort(faces,2);
faces_value=faces_temp(:,1)*pn^2+faces_temp(:,2)*pn+faces_temp(:,3);
[faces_value,index]=sort(faces_value);
surface_index=zeros(size(faces_value,1),1);
if faces_value(1)<faces_value(2)
    surface_index(1)=1;
end
if faces_value(size(faces_value,1))>faces_value(size(faces_value,1)-1)
    surface_index(size(faces_value,1))=1;
end
for i=2:size(faces_value,1)-1
    if faces_value(i)>faces_value(i-1) && faces_value(i+1)>faces_value(i)
        surface_index(i)=1;
    end
end
faces=faces(index,:);
t=faces(logical(surface_index),:);
p=x;

%%去除未使用的点
use=zeros(pn,1);
use(t(:))=1;
p=p(logical(use),:);
index=find(use);
use(index)=1:size(index,1);
t=use(t);