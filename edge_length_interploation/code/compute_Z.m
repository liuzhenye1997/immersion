%Z为绕Z轴旋转对应的旋转矩阵
function Z=compute_Z(theta)
Z=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0;0 0 1];


