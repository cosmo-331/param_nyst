function [Is] = param_nyst(As, ts, eps)
% Given an array of SPSD matrices As, compute an array of indices Is.
% Similar to AdaCUR
q = length(ts);
Is = cell(q,1);
[I, F] = RPCholesky(As(:,:,1), eps);
Is{1} = I;
%cnt = 0;
%h = 0;
for j = 2:q
    A = As(:,:,j);
    N = length(A);
    % Add new indices 
    [I, F, J] = RPCholesky_with_preindex(A, eps, Is{j-1});
    I_copy = I;
    
    % Remove redundant indices
    %k = length(I);
    if length(I) ~= length(Is{j-1})
        %h = h+1;
        for i = I_copy
            
            %R = chol(A(I,I));
            %C = A(:,I);
            %R_ = chol(A(setdiff(I,i), setdiff(I,i)));
            
            [I, J, error] = reduce_column_sketch(A,I,i,eps/N,J);

            %errors = [errors error];
            %true_error = norm(C/R, 'fro') - norm(A(:,setdiff(I,i))/R_, 'fro');
            %true_errors = [true_errors true_error];

        end
        %I = reduce_column(A,I,eps/N);
    end
    %if length(I)~= k
    %    cnt = cnt + 1;
    %end
    Is{j} = I;
end
%cnt
%h
%figure
%histogram(true_errors./errors);
%ylabel('frequency');
%label1 = xlabel('$\frac{Actual Error}{Estimated Error}$');
%set(label1, 'Interpreter','latex');