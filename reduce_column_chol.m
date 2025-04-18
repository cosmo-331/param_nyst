function [I] = reduce_column_chol(A,I,eps)
% Takes in SPSD A, U = pinv(A(I,I)), see if we need to remove i.
% Returns the updated U and I.
C = A(:,I);
C_squared = C'*C;
I_copy = I;
for i = I_copy
    r = length(I);
    j = find(I == i);
    R = chol(A(I,I));
    R_tilde = chol(A(setdiff(I,i), setdiff(I,i) ) ); 
    % This should take O(Nr) operations
    error = trace((R')\C_squared/(R)) - trace((R_tilde')\C_squared(setdiff(1:r,j), setdiff(1:r,j))/(R_tilde));
    if error < eps
        I = setdiff(I, i);
        C_squared = C_squared(setdiff(1:r,j), setdiff(1:r,j));
    end
end