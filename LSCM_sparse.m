%LSCM参数化
function new_points=LSCM_sparse(length,faces)
    xy=zeros(size(faces,1),6);
    for i=1:size(faces,1)
        theta=acos((length(i,1).^2+length(i,3).^2-length(i,2).^2)/(2*length(i,1)*length(i,3)));
        xy(i,3)=length(i,1);
        xy(i,4)=0;
        xy(i,5)=length(i,3)*cos(theta);
        xy(i,6)=length(i,3)*sin(theta);
    end

    point_number=max(faces(:));
    face_number=size(faces,1);

    %寻找边界固定点
    vertex_src=[faces(:,1);faces(:,2);faces(:,3)];
    vertex_dst=[faces(:,2);faces(:,3);faces(:,1)];
    edge=sort([vertex_src vertex_dst],2);
    [sedge,index]=sort(edge(:,1)*point_number+edge(:,2));
    edge=edge(index,:);
    % edge_temp=edge(1,:);
    % static;
    for i=1:size(sedge)-2
        if sedge(i+1)>sedge(i) && sedge(i+2)>sedge(i+1)
            static=edge(i+1,:);
            break;
        end
    end
    l=0;
    for i=1:face_number
        for j=1:3
            if faces(i,j)==static(1) && faces(i,mod(j,3)+1)==static(2)
                l=length(i,j);
            end
            if faces(i,j)==static(2) && faces(i,mod(j,3)+1)==static(1)
                l=length(i,j);
            end
        end
    end


    area=face_area(length);
    W_re=zeros(point_number,size(faces,1));
    W_im=zeros(point_number,size(faces,1));
    col=zeros(face_number,3);
    row=zeros(face_number,3);
    value_re=zeros(face_number,3);
    value_im=zeros(face_number,3);
    for i=1:size(faces,1)
        for j=1:3
            row(i,j)=faces(i,j);
            col(i,j)=i;
            value_re(i,j)=(-xy(i,2*(mod(j,3)+1)-1)+xy(i,2*(mod(j+1,3)+1)-1)); 
            value_im(i,j)=(-xy(i,2*(mod(j,3)+1))+xy(i,2*(mod(j+1,3)+1)));
        end
    end
    W_re=sparse(row(:),col(:),value_re(:),point_number,face_number);
    W_im=sparse(row(:),col(:),value_im(:),point_number,face_number);
    % l=compute_l(faces,points);
    % [min_l,path,path_length]=min_length(faces,l);
    % max_l=max(min_l(:));
    % index=find(min_l==max_l);
    % row=mod(index-1,point_number)+1;
    % col=(index-row)/point_number+1;
    % static=[row col];

    free=1:point_number;
    free(static)=[];
    M_re=(W_re'./sqrt(2*area));
    M_im=(W_im'./sqrt(2*area));
    p=2;
    Mf_re=M_re(:,free);
    Mf_im=M_im(:,free);
    Mp_re=M_re(:,static);
    Mp_im=M_im(:,static);
    A=[Mf_re -Mf_im;Mf_im Mf_re];
    u=[0;l;0;0];
    B=[Mp_re -Mp_im;Mp_im Mp_re];
    b=-B*u;
    AA=A'*A;
    x=(AA)\(A'*b);
    new_points=zeros(point_number,3);
    new_points(free,1)=x(1:point_number-p);
    new_points(free,2)=x(point_number-p+1:2*(point_number-p));
    new_points(static(1),1:2)=[0 0];
    new_points(static(2),1:2)=[l 0];
    d=A*x-b;
    sum(d(:).^2);
    % trimesh(faces, new_points(:,1), new_points(:,2), new_points(:,3), 'edgecolor', 'k'); axis off; axis equal; title('output');