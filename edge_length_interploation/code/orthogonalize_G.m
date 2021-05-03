%将G正交化
function G=orthogonalize_G(G)
for i=1:(size(G,1)/3)
    [U,S,V]=svd(G(3*i-2:3*i,:));
    if(det(U)*det(V)<0)
        V(3,:)=-V(3,:);
    end
    G(3*i-2:3*i,:)=U*V';
end
