%求得最优步长，用于四次能量
function alpha = best_step(points, faces, d, l_target)
% coefficient=zeros(4,1);
% for i=1:size(faces,1)
%     for j=1:3
%         dxi=d(faces(i,j),1);
%         dyi=d(faces(i,j),2);
%         dzi=d(faces(i,j),3);
%         dxj=d(faces(i,mod(j,3)+1),1);
%         dyj=d(faces(i,mod(j,3)+1),2);
%         dzj=d(faces(i,mod(j,3)+1),3);
%         xi=points(faces(i,j),1);
%         yi=points(faces(i,j),2);
%         zi=points(faces(i,j),3);
%         xj=points(faces(i,mod(j,3)+1),1);
%         yj=points(faces(i,mod(j,3)+1),2);
%         zj=points(faces(i,mod(j,3)+1),3);
%         l=l_target(i,j)^2;
%         coefficient(1,1)=coefficient(1,1)+2*((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2)*(2*(dxi - dxj)^2 + 2*(dyi - dyj)^2 + 2*(dzi - dzj)^2);
%         coefficient(2,1)=coefficient(2,1)+(2*((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2)*(2*(dxi - dxj)*(xi - xj) + 2*(dyi - dyj)*(yi - yj) + 2*(dzi - dzj)*(zi - zj)) + 2*(2*(dxi - dxj)*(xi - xj) + 2*(dyi - dyj)*(yi - yj) + 2*(dzi - dzj)*(zi - zj))*(2*(dxi - dxj)^2 + 2*(dyi - dyj)^2 + 2*(dzi - dzj)^2));
%         coefficient(3,1)=coefficient(3,1)+(2*(2*(dxi - dxj)*(xi - xj) + 2*(dyi - dyj)*(yi - yj) + 2*(dzi - dzj)*(zi - zj))^2 + 2*(2*(dxi - dxj)^2 + 2*(dyi - dyj)^2 + 2*(dzi - dzj)^2)*((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2));
%         coefficient(4,1)=coefficient(4,1)+2*(2*(dxi - dxj)*(xi - xj) + 2*(dyi - dyj)*(yi - yj) + 2*(dzi - dzj)*(zi - zj))*((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2);
%     end
% end
coefficient = best_step_c(points, faces, d, l_target);
root = -roots(coefficient);
root = real(root(imag(root) == 0));
alpha = min(root(root > 0));
if size(alpha, 1) == 0
    alpha = 0;
end