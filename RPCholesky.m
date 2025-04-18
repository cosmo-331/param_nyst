function [I, F, error, J] = RPCholesky(A, eps)
% Randomly Pivoted Cholesky. Keeps choosing index until trace norm of
% residual is less than eps.
N = length(A);
F = zeros(N,1);
d = diag(A);
I = [];
i = 1;
error = sum(d);
while error >= eps
    index = randsample(N,1,true,d);
    I = [I index];
    g = A(:,index);
    g = g - F(:, 1:i-1)*F(index, 1:i-1)';
    F(:,i) = g/sqrt(g(index));
    d = d - F(:,i).^2;
    error = error - norm(F(:,i))^2;
    d = max(d,0);
    % Sometimes the algorithm picks one index multiple times because of
    % floating point errors. Hence set d(i) = 0
    d(index) = 0;
    i = i+1;
end