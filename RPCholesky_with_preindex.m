function [I, F, J] = RPCholesky_with_preindex(A, eps, chosen_I)
% Randomly Pivoted Cholesky. Start from some given indices.
N = length(A);
R = chol(A(chosen_I, chosen_I));
F = A(:,chosen_I)/R;
d = diag(A);
I = chosen_I;
d = d - vecnorm(F,2,2).^2;
d(chosen_I) = 0;
i = length(chosen_I) + 1;
while sum(d) >= eps 
    d = max(d,0);
    index = randsample(N,1,true,d);
    I(i) = index;
    g = A(:,index);
    g = g - F(:, 1:i-1)*F(index, 1:i-1)';
    F(:,i) = g/sqrt(g(index));
    d = d - F(:,i).^2;
    % Sometimes the algorithm picks one index multiple times because of
    % floating point errors. Hence set d(i) = 0
    d(index) = 0;
    i = i+1;
end
J = zeros(i,1);
d = max(d,0);
for j = 1:i
    J(j) = randsample(1:N, 1, true, d);
    d(J(j)) = 0;
end