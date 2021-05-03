%a,b,c,d为二面角的四个点的坐标值
%theta则为需要旋转的角度
function theta=dihedral_angle(a,b,c,d)
v1=c-a;
v2=b-c;
n1=cross(v1,v2);
u1=d-c;
u2=b-d;
n2=cross(u1,u2);
n1=n1./norm(n1);
n2=n2./norm(n2);
if(dot(u1,n1)<0)
    theta=acos(dot(n1,n2)); 
else
    theta=-acos(dot(n1,n2)); 
end
