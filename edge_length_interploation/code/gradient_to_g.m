%gradient为矩阵形式，gradient(i,j)为点i和点j连线对应的二面角的梯度，方便检索
%g为向量形式，即将梯度排成一列，方便运算
function g=gradient_to_g(faces,gradient,A_index)
face_number=size(faces,1);
g=zeros(3*face_number/2,1);
[row,col,val]=find(A_index);
index=sub2ind([size(A_index,1) size(A_index,1)],row,col);
g(val)=gradient(index)*2;



