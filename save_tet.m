function save_tet(x,t,name)
fid=fopen(name,'w');
fprintf(fid,'MeshVersionFormatted 1\n');
fprintf(fid,'Dimension\n3\n\n');
x=x';
fprintf(fid,'Vertices\n%d\n',size(x,2));
fprintf(fid,'%.20f %.20f %.20f\r\n',x(:)); 
fprintf(fid,'\nTriangles\n0\n\n');
t=t';
fprintf(fid,'Tetrahedra\n%d\n',size(t,2));
fprintf(fid,'%d %d %d %d\r\n',t(:)); 
fclose(fid);