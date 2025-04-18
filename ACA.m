function [i] = ACA(A,I)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
R = chol(A(I,I));
F = A(:,I)/R;
d = diag(A);
for i = 1:length(I)
    d = d - F(:,i).^2;
end
[aaaa, i] = max(d);
end