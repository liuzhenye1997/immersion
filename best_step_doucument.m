syms alpha dxi dxj  dyi dyj dzi dzj xi xj yi yj zi zj l
%��ʽΪ��x+alpha*d��������E�ڱ�(i,j)�ϵ�ֵ
Eij = ((xi - xj + alpha * (dxi - dxj))^2 + (yi - yj + alpha * (dyi - dyj))^2 + (zi - zj + alpha * (dzi - dzj))^2 - l)^2;
%��alpha�󵼺󻯼�
dEij = diff(Eij, alpha);
collect(dEij, alpha)
%�����õ���ʽ����һ������alpha�����κ���
2 * ((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2) * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2) * alpha^3 + (2 * ((dxi - dxj)^2 + (dyi - dyj)^2 + (dzi - dzj)^2) * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) + 2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2)) * alpha^2 + (2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj))^2 + 2 * (2 * (dxi - dxj)^2 + 2 * (dyi - dyj)^2 + 2 * (dzi - dzj)^2) * ((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2)) * alpha + 2 * (2 * (dxi - dxj) * (xi - xj) + 2 * (dyi - dyj) * (yi - yj) + 2 * (dzi - dzj) * (zi - zj)) * ((xi - xj)^2 - l + (yi - yj)^2 + (zi - zj)^2)
%�����еı��ϵ�ϵ�����,�õ�dE=a*alpha^3+b*alpha^2+c*alpha+d��
%��dE=0ʱ��E�ﵽ���ֵ�������dE==0���һԪ���η��̼��ɵõ���Ѳ���alpha