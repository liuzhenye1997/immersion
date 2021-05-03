%求点的度
function degree = points_degree(faces)
face_number=size(faces,1);
point_number = max(faces(:));
A = zeros(point_number, point_number);
for i = 1:face_number
    for j = 1:3
        A(faces(i, j), faces(i, mod(j, 3) + 1)) = 1;
        A(faces(i, j), faces(i, mod(j - 2, 3) + 1)) = 1;
    end
end
degree = sum(A);
