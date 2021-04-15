%计算能量Ef所用矩阵
function A=conpute_Ef_matrix(faces,R)
    face_number=size(faces,1);
    ad_face=sparse(faces,[faces(size(faces,1)+1:3*size(faces,1)),faces(1:size(faces,1))]',[1:face_number 1:face_number 1:face_number]);
 
    col=zeros(face_number,3,3,3);
    row=zeros(face_number,3,3,3);
    val=zeros(face_number,3,3,3);
    
    index1=[1 2 3;1 2 3;1 2 3];
    index2=[1 1 1;2 2 2;3 3 3];
    for i=1:face_number
        for j=1:3
            f1=ad_face(faces(i,j),faces(i,mod(j+1,3)+1));
            col(i,j,:,:)=9*(i-1)+3*(j-1)+ index1;
            row(i,j,:,:)=3*(f1-1)+index2;
        end
    end
    %A1是能量中fi*Rij部分的系数
    A1=sparse(row(:),col(:),R(:),3*face_number,9*face_number);

    I=eye(3);
    for i=1:face_number
        for j=1:3
            f2=ad_face(faces(i,mod(j+1,3)+1),faces(i,j));
            row(i,j,:,:)=3*(f2-1)+index2;
            val(i,j,:,:)=I;
        end
    end
    %A2是能量中fj部分的系数
    A2=sparse(row(:),col(:),val(:),3*face_number,9*face_number);
    A=A1-A2;