function [I, J, error] = reduce_column_sketch(A,I,i,eps,J)

j = find(I==i);
K = union(J, i);
R = chol(A(setdiff(I,i),setdiff(I,i)));
F = A(K,setdiff(I,i))/R;
error = norm(F, 'fro')^2;
if error < eps
    I = setdiff(I,i);
end