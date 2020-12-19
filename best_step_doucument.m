syms alpha dxi dxj  dyi dyj dzi dzj xi xj yi yj zi zj l
%下式为在x+alpha*d处的能量E在边(i,j)上的值
Eij = ((xi - xj + alpha * (dxi - dxj))^2 + (yi - yj + alpha * (dyi - dyj))^2 + (zi - zj + alpha * (dzi - dzj))^2 - l)^2;
%对alpha求导后化简
dEij = diff(Eij, alpha);
collect(dEij, alpha)
%化简后得到下式，是一个关于alpha的三次函数
2 * ((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2) * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2) * alpha^3 + (2 * ((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2) * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) + 2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2)) * alpha^2 + (2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj))^2 + 2 * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2) * ((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2)) * alpha + 2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) * ((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2)
%将所有的边上的系数相加,得到dE=a*alpha^3+b*alpha^2+c*alpha+d。
%当dE=0时，E达到最大值，故求解dE==0这个一元三次方程即可得到最佳步长alpha