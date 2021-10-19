function [X, out] = globalCGLS(A, B, param)
% [X, out] = globalCGLS(A, B, param)
% applies the global CGLS method to a least squares problem with 
% multiple-right hand sides min_X || B - A X||_F
% 
% inputs
% A, B: complex matrices
% param.maxit: maximum number of iterations
% param.tol: stopping criterion for the iterations
% 
% outputs
% X: approximate solution
% out.iter: required number of iterations
% out.resnrm: relative residual norm
% 
% reference
% Xiao-Guan Lv, Ting-Zhu Huang, Le Jiang, Jun Liu,
% Two global iterative methods for ill-posed problems from image restoration
% Journal of Computational Analysis and Applications
% Volume 18, Number 2, pp. 219-237, 2015.

maxit = param.maxit;
tol = param.tol;

resnrm = zeros(maxit, 1);
[~, n] = size(A);
l = size(B, 2);

X = zeros(n, l);

AT = A';
R = B;
S = AT * R;
P = S;

nrmATR0 = norm(S, 'fro');
FnrmS_old = nrmATR0;

% main loop
for k = 1:maxit
    Q = A * P;
    
    FnrmQ = norm(Q, 'fro');
    
    alp = FnrmS_old*FnrmS_old / (FnrmQ * FnrmQ);
    
    X = X + alp*P;
    
    R = R - alp*Q;
    
    S = AT * R;
    
    FnrmS_new = norm(S, 'fro');

    resnrm(k) = FnrmS_new / nrmATR0;
    
%   stopping rule
    if resnrm(k) <= tol
        out.iter = k;
        out.resnrm = resnrm(1:k);
        return
    end
    
    beta = FnrmS_new*FnrmS_new / (FnrmS_old*FnrmS_old);
    FnrmS_old = FnrmS_new;
    
    P = S + beta*P;
end

fprintf('globalCGLS failded to converge\n');
out.resnrm = resnrm(1:k-1);
out.iter = k-1;