function I = param_ACA(As,phis, eps)
% param_ACA in Kressner's paper
% phis is a matrix where 
% phis(i,j) = phi_i(theta_j)
I = [];
s = size(As,3);
N = size(As,1);
traces = zeros(s,1);
m = size(phis,2);
for j = 1:s
    traces(j) = trace(As(:,:,j));
end
while 0<1
    theta_star = 1;
    r = length(I);
    res_max = 0;
    temp = zeros(N, s*r);
    
    if r > 0
        for j=1:s
            [Q_I, R_I] = insertcolumn(Q_I,R_I,As(:,I(r),j),j*r);
        end
    else
        for j = 1:s
            temp(:,j*r-r+1: j*r) = As(:,I,j);
        end
        [Q_I, R_I] = qr(temp, "econ");
    end
    for i=1:m
        A_I = zeros(r);
        for j=1:s
            A_I = A_I + phis(j,i)*As(I,I,j);
        end
        %min(eig(A_I))
        R_A = chol(A_I);
        res = 0;
        for j=1:s
            res = res + phis(j,i)*traces(j);
        end
        phi_theta = zeros(r*s, r);
        for j = 1:s
            phi_theta(j*r - r +1:j*r,:) = eye(r)*phis(j,i);
        end
        res = res - norm(R_I*phi_theta/R_A, 'fro')^2;
        if res >= res_max
            res_max = res;
            theta_star = i;
        end
    end
    if res_max <= eps
        break
    end
    i_new = ACA(As(:,:,theta_star), I);
    I = [I i_new];
end
end