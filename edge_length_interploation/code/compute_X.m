%X为绕X轴旋转对应的旋转矩阵
function X=compute_X(theta)
X=[1 0 0;0 cos(theta) -sin(theta);0 sin(theta) cos(theta)];