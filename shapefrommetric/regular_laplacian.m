%����laplacian����d0ÿһ�д���һ���ߣ�-1��������㣬1��������㡣
%��d0'*d0��Ϊһ���laplacian��������ʹ����point_weight�����˼�Ȩ��
function [L,M,d0,star1]=regular_laplacian(face_weight,point_weight,face_number,point_number,vertex_src,vertex_dst)
vertex_src=vertex_src';
vertex_dst=vertex_dst';
vertex_src=vertex_src(:);
vertex_dst=vertex_dst(:);
d0_row=[1:3*face_number,1:3*face_number];
d0_col=[vertex_src,vertex_dst];
d0_val=[-ones(face_number*3,1),ones(face_number*3,1)];
d0=sparse(d0_row,d0_col,d0_val,face_number*3,point_number);

star0_right=sparse(1:length(point_weight(:)),1:length(point_weight(:)),point_weight);
star1_right=sparse(1:3*face_number,1:3*face_number,face_weight);
L=-d0'*(star1_right*d0);
M=star0_right;
star1=star1_right;