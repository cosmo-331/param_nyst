function [I] = reduce_column(A,I,eps)
% Takes in SPSD A, U = pinv(A(I,I)), see if we need to remove i.
% Returns the updated U and I.
I = setdiff(I,-1);
C = A(:,I);
U = inv(A(I,I));
I_copy = I;
for i = I_copy
    r = length(I);
    X = zeros(r,2);
    j = find(I == i);
    X(j,1) = 1;
    X(:,2) = A(I,i);
    X(j,2) = 0;
    W = X';
    W([1,2],:) = W([2,1],:);
    % This should take O(Nr) operations
    error = norm(A(:,i))^2/A(i,i) - trace(((eye(2) - W*U*X)\W*U*C')*(C*(U*X)));
    if error < eps
        I = setdiff(I, i);
        U = U + (U*X/(eye(2) - W*U*X))*W*U;
        U = U(setdiff(1:r, j), setdiff(1:r,j));
        C = C(:,setdiff(1:r, j));
    end
end