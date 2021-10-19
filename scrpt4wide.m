% generate A and B of matrix pencil zB-A, m < n
clear; close all; rng(0, 'twister');

% uncomment one of the following cases
m = 10; n = 100; nu = 0;
% m = 100; n = 1000; nu = 0;
% m = 1000; n = 10000; nu = 0;
% m = 10000; n = 100000; nu = 0;

center = 1 + 1i; % center of the circle region

% set the radius of the circle region for each (m, n)
switch n
    case 100
        if m == 10; radius = 1; end
    case 1000
        if m == 100; radius = 0.3; end
    case 10000
        if m == 1000; radius = 0.1; end
    case 100000
        if m == 10000; radius = 0.05; end
    otherwise
        warning('n is wrong')
end

fprintf('m = %d, n = %d, nu = %d, radius = %.2f\n', 3*m, n, nu, radius);

%% finite eigenvalues
lmd = randn(m, 1) + 1i*randn(m, 1); % exact eigenvalue
t = length(lmd(abs(lmd-center) < radius));
fprintf('# of target eigenvalues: %d\n', t);

figure(1) 
fig0 = plot(lmd, 'go');
hold on

if m == 10000 && n == 100000
    load A30000x100000
    load B30000x100000
    load lmd30000x100000
    return
end

A1 = sparse(1:m, 1:m, lmd, m, m);
B1 = sparse(1:m, 1:m, ones(m, 1), m, m);

%% infinite eigenvalues
r_ind = sort(randperm(m-1, m/2));
A2 = speye(m, m);
B2 = sparse(r_ind+1, r_ind, 1, m, m);

A = [[A1, sparse(m, m);
     sparse(m, m), A2], sparse(2*m, n-2*m);
     sparse(m, n)];
B = [[B1, sparse(m, m);
     sparse(m, m), B2], sparse(2*m, n-2*m);
     sparse(m, n)];

clear A1 B1 A2 B2

if n > 1000
    tic
    while min(nnz(A), nnz(B)) < .95*m*n*0.001
        [A, B] = myrjr2(A, B, 2);
    end
    toc
    return
end

A = full(A);
B = full(B);

R = randn(n)+1i*randn(n);
A = A * R;
B = B * R;

R = randn(3*m)+1i*randn(3*m);
A = R * A;
B = R * B;