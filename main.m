%% Generate Numerical Examples
N = 1000;
T = 101;
A = randn(N, N);
A_skew = tril(A,-1) - triu(A', 1);
D = diag(1./2.^(1:N));
ts = linspace(0,1,T);
As = zeros(N,N,T);
for i = 1:T
    U = expm(A_skew*ts(i));
    As(:,:,i) = U*exp(ts(i))*D*U';
end

%% Run the algorithm and compute the errors
figure
tic
Is = param_nyst(As,ts,1e-6);
toc
ranks = zeros(T,3);
true_ranks = zeros(T,3);
errors = zeros(T,1);
for i = 1:T
    I = Is{i};
    A = As(:,:,i);
    ranks(i,1) = length(I);
    true_ranks(i,1) = rank(A, 1e-9);
    R = chol(A(I,I));
    %A_hat = (A(:,I)/(R))*(A(:,I)/(R))';
    errors(i) = trace(A) - norm(A(:,I)/(R),'fro')^2;
end
semilogy(ts, errors, '-b.');

hold on
tic
Is = param_nyst(As,ts,1e-8);
toc
errors = zeros(T,1);
for i = 1:T
    I = Is{i};
    A = As(:,:,i);
    ranks(i,2) = length(I);
    true_ranks(i,2) = rank(A, 1e-11);
    R = chol(A(I,I));
    A_hat = (A(:,I)/(R))*(A(:,I)/(R))';
    errors(i) = trace(A) - norm(A(:,I)/(R),'fro')^2;
end
semilogy(ts, errors, '-r.');

tic
Is = param_nyst(As,ts,1e-10);
toc
errors = zeros(T,1);
for i = 1:T
    I = Is{i};
    A = As(:,:,i);
    ranks(i,3) = length(I);
    true_ranks(i,3) = rank(A, 1e-13);
    R = chol(A(I,I));
    A_hat = (A(:,I)/(R))*(A(:,I)/(R))';
    errors(i) = trace(A) - norm(A(:,I)/(R),'fro')^2;
end

semilogy(ts, errors, '-g.');

leg1 = legend('$\varepsilon = 10^{-6}$','$\varepsilon = 10^{-8}$','$\varepsilon = 10^{-10}$', Location='southeast');
set(leg1,'Interpreter','latex');
xlabel('Parameter t')
ylabel('Trace norm error')
grid on
hold off

figure
plot(ts, ranks(:,1), '-b.');
hold on
plot(ts, ranks(:,2), '-r.');
plot(ts, ranks(:,3), '-g.');
xlabel('Parameter t')
ylabel('Rank')
grid on
hold off
%%
plot(ts, true_ranks(:,1), '-.', 'Color','black');
plot(ts, true_ranks(:,2), '-.', 'Color','black');
plot(ts, true_ranks(:,3), '-.', 'Color','black');
leg1 = legend('$\varepsilon = 10^{-6}$','$\varepsilon = 10^{-8}$','$\varepsilon = 10^{-10}$', '$\mathtt{rank}_{\varepsilon/N}(\mathbf{A}(t))$', Location='southeast');
set(leg1,'Interpreter','latex');
hold off

%% RPCholesky runtime
tic
for i=1:T
    I = RPCholesky(As(:,:,i),1e-10);
end
toc

%% Compute index only once
figure
%[I, F] = RPCholesky(As(:,:,1), 1e-9);
%r = length(I);
errors = zeros(T,1);
for i = 1:T
    A = As(:,:,i);
    ranks(i) = length(I);
    R = chol(A(I,I));
    A_hat = (A(:,I)/(R))*(A(:,I)/(R))';
    errors(i) = trace(A) - trace(A_hat);
    errors(i) = errors(i)/exp(ts(i));
end
bounds = zeros(T,1);
for i = 1:T
    A = As(:,:,i);
    [Q_C,R] = qr(A(:,I), 'econ');
    V_r = expm(A_skew*ts(i));
    V_r = V_r(:,1:r);
    bounds(i) = norm(pinv(Q_C(I,:)))*norm(pinv(V_r(I,:)))/2^r/exp(ts(i));
end
semilogy(ts(2:T), errors(2:T), '-b.')
hold on
semilogy(ts(2:T), bounds(2:T), '-g.')
hold off
xlabel('Parameter t')
leg1 = legend('$\displaystyle \frac{||\mathbf{A} - \mathbf{C}\mathbf{U}^{-1}\mathbf{C}^T||_*}{||\mathbf{A}||_*}$','$\displaystyle \frac{||\mathbf{Q}_{\mathbf{C}}(I,:)^\dagger||_2||\mathbf{V}_{r}(I,:)^{-1}||_2||\mathbf{A} - [\mathbf{A}]_r||_*}{||\mathbf{A}||_*}$', Location='southeast');
set(leg1,'Interpreter','latex');
grid on

%% Rank-adaptivity Test 
% 
N = 1000;
T = 10;

ts = 1:T;
As = zeros(N,N,T);
temp = randn(N,N);
[Q,qwe] = qr(temp, 'econ');
D = diag(1./2.^(1:N));
As(:,:,1) = Q*D*Q';
for i = 2:T
    G = randn(N,2)/sqrt(N)/100;
    As(:,:,i) = G*G' + As(:,:,i-1);
end

phis = eye(T);
tic
I = param_ACA(As,phis,1e-8);
toc
errors = zeros(T,3);
ranks = zeros(T,3);
r = size(I);
for i = 1:T
    A = As(:,:,i);
    ranks(i,1) = length(I);
    %true_ranks(i,1) = rank(A, 1e-9);
    R = chol(A(I,I));
    errors(i,1) = trace(A) - norm(A(:,I)/R, 'fro')^2;
end
figure
semilogy(ts, errors(:,1), '-b.');

Is = param_nyst(As,ts,1e-8);
for i = 1:T
    I = Is{i};
    A = As(:,:,i);
    ranks(i,2) = length(I);
    %true_ranks(i,1) = rank(A, 1e-9);
    R = chol(A(I,I));
    errors(i,2) = trace(A) - norm(A(:,I)/R, 'fro')^2;
end

hold on
semilogy(ts, errors(:,2), '-r.');

for i = 1:T
    A = As(:,:,i);
    I = RPCholesky(A,1e-8);
    ranks(i,3) = length(I);
    %true_ranks(i,1) = rank(A, 1e-9);
    R = chol(A(I,I));
    errors(i,3) = trace(A) - norm(A(:,I)/R, 'fro')^2;
end
semilogy(ts, errors(:,3), '-m.');
legend('Param\_ACA','Param\_Nyst','RPCholesky',location='southeast');
xlabel('Parameter t');
ylabel('Trace norm error');
hold off

figure
plot(ts, ranks(:,1), '-b.');
hold on
plot(ts, ranks(:,2), '-r.');
plot(ts, ranks(:,3), '-m.');
legend('Param\_ACA','Param\_Nyst','RPCholesky',location='southeast');
xlabel('Parameter t');
ylabel('Rank');
hold off

%% Speed Test (adversarial)
N = 5000;
temp = randn(N,400);
[U, temp] = qr(temp, 'econ');
D  = diag(logspace(-16,0,400));
evalin('base', 'clear temp');
T = 40;
ts = 1:T;
As = zeros(N,N,T);
As(:,:,1) = U*D*U';
for i=2:T
    G = randn(5000,4)/sqrt(5000)/100;
    As(:,:,i) = As(:,:,i-1) + G*G';
end
%% Our Algorithm
for p=1:10
tic
Is = param_nyst(As,ts,1e-6);
toc
end
%% RPCholesky
for p=1:10
tic
for i=1:T
    A = As(:,:,i);
    I = RPCholesky(A,1e-6);
end
toc
end
%% param_ACA (too slow)
phis = eye(T);
tic
I = param_ACA(As(:,:,1:20), phis(1:20,1:20),1e-6);
toc