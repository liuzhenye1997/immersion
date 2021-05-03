%A为矩阵形式，A(i,j)为点i和点j连线对应的二面角，方便检索
%a为向量形式，即将二面角排成一列，方便运算
function A=atoA(faces,A_index,a)
[row,col,val]=find(A_index);
A=sparse(row,col,a(val));