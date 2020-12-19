%求能量的hessian矩阵
function H = hessian_matrix_sparse(points, faces, l_target)
% point_number=size(points,1);
% col=zeros(size(faces,1),3,36);
% row=zeros(size(faces,1),3,36);
% value=zeros(size(faces,1),3,36);
% for i=1:size(faces,1)
%     for j=1:3
%         xi=points(faces(i,j),1);
%         yi=points(faces(i,j),2);
%         zi=points(faces(i,j),3);
%         xj=points(faces(i,mod(j,3)+1),1);
%         yj=points(faces(i,mod(j,3)+1),2);
%         zj=points(faces(i,mod(j,3)+1),3);
%         xi_index=faces(i,j);
%         yi_index=faces(i,j)+point_number;
%         zi_index=faces(i,j)+2*point_number;
%         xj_index=faces(i,mod(j,3)+1);
%         yj_index=faces(i,mod(j,3)+1)+point_number;
%         zj_index=faces(i,mod(j,3)+1)+2*point_number;
%         l=l_target(i,j);
% %         H(xi_index,xi_index)=H(xi_index,xi_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*xi - 2*xj)^2 - 4*l^2;
% %         H(xj_index,xj_index)=H(xj_index,xj_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*xi - 2*xj)^2 - 4*l^2;
% %         H(yi_index,yi_index)=H(yi_index,yi_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*yi - 2*yj)^2 - 4*l^2;
% %         H(yj_index,yj_index)=H(yj_index,yj_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*yi - 2*yj)^2 - 4*l^2;
% %         H(zi_index,zi_index)=H(zi_index,zi_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*zi - 2*zj)^2 - 4*l^2;
% %         H(zj_index,zj_index)=H(zj_index,zj_index)+4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*zi - 2*zj)^2 - 4*l^2;
%         row(i,j,1)=xi_index;col(i,j,1)=xi_index;value(i,j,1)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*xi - 2*xj)^2 - 4*l^2;
%         row(i,j,2)=xj_index;col(i,j,2)=xj_index;value(i,j,2)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*xi - 2*xj)^2 - 4*l^2;
%         row(i,j,3)=yi_index;col(i,j,3)=yi_index;value(i,j,3)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*yi - 2*yj)^2 - 4*l^2;
%         row(i,j,4)=yj_index;col(i,j,4)=yj_index;value(i,j,4)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*yi - 2*yj)^2 - 4*l^2;
%         row(i,j,5)=zi_index;col(i,j,5)=zi_index;value(i,j,5)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*zi - 2*zj)^2 - 4*l^2;
%         row(i,j,6)=zj_index;col(i,j,6)=zj_index;value(i,j,6)=4*(xi - xj)^2 + 4*(yi - yj)^2 + 4*(zi - zj)^2 + 2*(2*zi - 2*zj)^2 - 4*l^2;
%
%
%
% %         H(xi_index,yi_index)=H(xi_index,yi_index)+2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(yi_index,xi_index)=H(yi_index,xi_index)+2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(zi_index,xi_index)=H(zi_index,xi_index)+2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(xi_index,zi_index)=H(xi_index,zi_index)+2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(yi_index,zi_index)=H(yi_index,zi_index)+2*(2*yi - 2*yj)*(2*zi - 2*zj);
% %         H(zi_index,yi_index)=H(zi_index,yi_index)+2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,7)=xi_index;col(i,j,7)=yi_index;value(i,j,7)=2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,8)=yi_index;col(i,j,8)=xi_index;value(i,j,8)=2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,9)=zi_index;col(i,j,9)=xi_index;value(i,j,9)=2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,10)=xi_index;col(i,j,10)=zi_index;value(i,j,10)=2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,11)=yi_index;col(i,j,11)=zi_index;value(i,j,11)=2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,12)=zi_index;col(i,j,12)=yi_index;value(i,j,12)=2*(2*yi - 2*yj)*(2*zi - 2*zj);
%
%
% %         H(xi_index,xj_index)=H(xi_index,xj_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*xi - 2*xj)^2 - 4*(xi - xj)^2;
% %         H(xj_index,xi_index)=H(xj_index,xi_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*xi - 2*xj)^2 - 4*(xi - xj)^2;
% %         H(yi_index,yj_index)=H(yi_index,yj_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*yi - 2*yj)^2 - 4*(xi - xj)^2;
% %         H(yj_index,yi_index)=H(yj_index,yi_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*yi - 2*yj)^2 - 4*(xi - xj)^2;
% %         H(zi_index,zj_index)=H(zi_index,zj_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*zi - 2*zj)^2 - 4*(xi - xj)^2;
% %         H(zj_index,zi_index)=H(zj_index,zi_index)+4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*zi - 2*zj)^2 - 4*(xi - xj)^2;
%         row(i,j,13)=xi_index;col(i,j,13)=xj_index;value(i,j,13)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*xi - 2*xj)^2 - 4*(xi - xj)^2;
%         row(i,j,14)=xj_index;col(i,j,14)=xi_index;value(i,j,14)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*xi - 2*xj)^2 - 4*(xi - xj)^2;
%         row(i,j,15)=yi_index;col(i,j,15)=yj_index;value(i,j,15)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*yi - 2*yj)^2 - 4*(xi - xj)^2;
%         row(i,j,16)=yj_index;col(i,j,16)=yi_index;value(i,j,16)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*yi - 2*yj)^2 - 4*(xi - xj)^2;
%         row(i,j,17)=zi_index;col(i,j,17)=zj_index;value(i,j,17)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*zi - 2*zj)^2 - 4*(xi - xj)^2;
%         row(i,j,18)=zj_index;col(i,j,18)=zi_index;value(i,j,18)=4*l^2 - 4*(yi - yj)^2 - 4*(zi - zj)^2 - 2*(2*zi - 2*zj)^2 - 4*(xi - xj)^2;
%
% %         H(xi_index,yj_index)=H(xi_index,yj_index)-2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(yj_index,xi_index)=H(yj_index,xi_index)-2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(zj_index,xi_index)=H(zj_index,xi_index)-2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(xi_index,zj_index)=H(xi_index,zj_index)-2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,19)=xi_index;col(i,j,19)=yj_index;value(i,j,19)=-2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,20)=yj_index;col(i,j,20)=xi_index;value(i,j,20)=-2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,21)=zj_index;col(i,j,21)=xi_index;value(i,j,21)=-2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,22)=xi_index;col(i,j,22)=zj_index;value(i,j,22)=-2*(2*xi - 2*xj)*(2*zi - 2*zj);
%
% %         H(yi_index,xj_index)=H(yi_index,xj_index)-2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(xj_index,yi_index)=H(xj_index,yi_index)-2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(yi_index,zj_index)=H(yi_index,zj_index)-2*(2*yi - 2*yj)*(2*zi - 2*zj);
% %         H(zj_index,yi_index)=H(zj_index,yi_index)-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,23)=yi_index;col(i,j,23)=xj_index;value(i,j,23)=-2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,24)=xj_index;col(i,j,24)=yi_index;value(i,j,24)=-2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,25)=yi_index;col(i,j,25)=zj_index;value(i,j,25)=-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,26)=zj_index;col(i,j,26)=yi_index;value(i,j,26)=-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%
%
% %         H(xj_index,zi_index)=H(xj_index,zi_index)-2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(zi_index,xj_index)=H(zi_index,xj_index)-2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(yj_index,zi_index)=H(yj_index,zi_index)-2*(2*yi - 2*yj)*(2*zi - 2*zj);
% %         H(zi_index,yj_index)=H(zi_index,yj_index)-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,27)=xj_index;col(i,j,27)=zi_index;value(i,j,27)=-2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,28)=zi_index;col(i,j,28)=xj_index;value(i,j,28)=-2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,29)=yj_index;col(i,j,29)=zi_index;value(i,j,29)=-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,30)=zi_index;col(i,j,30)=yj_index;value(i,j,30)=-2*(2*yi - 2*yj)*(2*zi - 2*zj);
%
% %         H(xj_index,yj_index)=H(xj_index,yj_index)+2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(yj_index,xj_index)=H(yj_index,xj_index)+2*(2*xi - 2*xj)*(2*yi - 2*yj);
% %         H(xj_index,zj_index)=H(xj_index,zj_index)+2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(zj_index,xj_index)=H(zj_index,xj_index)+2*(2*xi - 2*xj)*(2*zi - 2*zj);
% %         H(zj_index,yj_index)=H(zj_index,yj_index)+2*(2*yi - 2*yj)*(2*zi - 2*zj);
% %         H(yj_index,zj_index)=H(yj_index,zj_index)+2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,31)=xj_index;col(i,j,31)=yj_index;value(i,j,31)=2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,32)=yj_index;col(i,j,32)=xj_index;value(i,j,32)=2*(2*xi - 2*xj)*(2*yi - 2*yj);
%         row(i,j,33)=xj_index;col(i,j,33)=zj_index;value(i,j,33)=2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,34)=zj_index;col(i,j,34)=xj_index;value(i,j,34)=2*(2*xi - 2*xj)*(2*zi - 2*zj);
%         row(i,j,35)=zj_index;col(i,j,35)=yj_index;value(i,j,35)=2*(2*yi - 2*yj)*(2*zi - 2*zj);
%         row(i,j,36)=yj_index;col(i,j,36)=zj_index;value(i,j,36)=2*(2*yi - 2*yj)*(2*zi - 2*zj);
%     end
% end
% [row,col,value]=hessian_matrix_triplet(points,faces,l_target);
% H=sparse(row(:),col(:),value(:));
[~, ~, ~, H] = hessian_matrix_sparse_c(points, faces, l_target);
