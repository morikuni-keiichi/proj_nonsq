% driver

[m, n] = size(A);

% set the values of parameters
if (m == 100 && n == 30) || (m == 30 && n == 100) || ...
   (m == 1000 && n == 300) || (m == 300 && n == 1000)
    param.L = 4;
    param.M = 2;
elseif (m == 10000 && n == 3000) || (m ==  3000 && n == 10000) || ...
       (m == 10000 && n == 2001) || (m == 10000 && n ==  2010) || ...
       (m == 10000 && n == 2100) || (m == 10000 && n ==  3000) || ...
       (m == 100000 && n == 30000) || (m == 30000 && n == 100000)
    param.L = 8;
    param.M = 4;
    param.tol = 1e-14;   
    param.maxit = min(size(A));
else
    error('m or n is wrong')
end

param.N = 48;
param.center = center;
param.radius = radius;

% Plot the target region
z = param.center + param.radius*exp(2*pi*1i*((0:param.N-1)+0.5)/param.N);
plot(z, 'm');

%% Apply the method to matrix pencil zB-A
fprintf('\n-- running --\n');
tic; [V2, lmd2] = proj_nonsq(A, B, param); toc

%% Plot the eigenvalues 
ind2 = abs(lmd2-param.center) < param.radius; 
fig2 = plot(lmd2(ind2), 'rx');
legend('Exact eigenvalue', 'Target region', 'Computed eigenvalue')

% Compute the relative resitual norm
nrmA = norm(A, 'fro');
nrmB = norm(B, 'fro');
resvec2 = zeros(size(V2, 2), 1);
for k = 1:size(V2, 2)
    resvec2(k) = norm((A - lmd2(k)*B)*V2(:, k)) / (nrmA+abs(lmd2(k))*nrmB);
end

% Compute the relative error
ind = find(abs(lmd-param.center) < param.radius);
for i = 1:length(ind)
    j = ind(i);
    [minval(i), ind3(i)] = min(abs(lmd2 - lmd(j)));
end

fprintf('relative residual norm %.2e for lmd = %.8e + %.8e\n', [resvec2(ind3)'; real(lmd2(ind3))'; imag(lmd2(ind3))']);
fprintf('relative error %.2e for lmd = %.8e + %.8e\n', [minval ./ abs(lmd(ind))'; real(lmd2(ind3))'; imag(lmd2(ind3))']);