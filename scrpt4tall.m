% generate A and B of matrix pencil zB-A, m > n
clear; close all; rng(0, 'twister');

% uncomment one of the following cases
m = 100; n = 10; nu = 0;
% m = 1000; n = 100; nu = 0;
% m = 10000; n = 1000; nu = 0;
% m = 100000; n = 10000; nu = 0;
% m = 10000; n = 1000; nu = 1;
% m = 10000; n = 1000; nu = 10;
% m = 10000; n = 1000; nu = 100;
% m = 10000; n = 1000; nu = 1000;

center = 1 + 1i; % center of the circle region

% set the radius of the circle region for each (m, n)
switch m
    case 100
        if n == 10; radius = 1; end 
    case 1000
        if n == 100; radius = 0.3; end
    case 10000
        if n == 1000; radius = 0.1; end
    case 100000
        if n == 10000; radius = 0.05; end
    otherwise
        warning('n is wrong')
end

fprintf('m = %d, n = %d, nu = %d, radius = %.2f\n', m, 3*n, nu, radius);

% finite eigenvalues
lmd = randn(n, 1) + 1i*randn(n, 1); % exact eigenvalue
t = length(lmd(abs(lmd-center) < radius));
fprintf('# of target eigenvalues: %d\n', t);

figure(1) 
fig0 = plot(lmd, 'go');
hold on
xlabel('Re', 'Interpreter', 'latex')
ylabel('Im', 'Interpreter', 'latex')


if m == 100000 && n == 10000
    load A100000x30000
    load B100000x30000
    load lmd100000x30000
    return
end

A1 = sparse(1:n, 1:n, lmd, n, n);
B1 = sparse(1:n, 1:n, ones(n, 1), n, n);

% infinite eigenvalues
r_ind = sort(randperm(n-1, n/2));
A2 = speye(n, n);
B2 = sparse(r_ind, r_ind+1, 1, n, n);

if nu == 0
    A = [A1, sparse(n, 2*n);
        sparse(n, n), A2, sparse(n, n);
        sparse(m-2*n, 3*n)];
    B = [B1, sparse(n, 2*n);
        sparse(n, n), B2, sparse(n, n);
        sparse(m-2*n, 3*n)];
elseif nu > 0
    A = [A1, sparse(n, n+nu);
        sparse(n, n), A2, sparse(n, nu);
        sparse(m-2*n, 2*n), sparse(2:nu+1, 1:nu, 1, m-2*n, nu)];
    B = [B1, sparse(n, n+nu);
        sparse(n, n), B2, sparse(n, nu);
        sparse(m-2*n, 2*n), sparse(1:nu, 1:nu, 1, m-2*n, nu)];
end
 
clear A1 B1 A2 B2

% generate sparse matrices
if m > 1000
    while min(nnz(A), nnz(B)) < .95*m*n*0.001
        [A, B] = myrjr2(A, B, 2);
    end    
    return
end

% generate dense matrices
A = full(A);
B = full(B);

R = randn(3*n+nu)+1i*randn(3*n+nu);
A = A * R;
B = B * R;

R = randn(m)+1i*randn(m);
A = R * A;
B = R * B;
