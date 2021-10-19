function [X, lmd] = proj_nonsq(A, B, param)
% [V, lmd] = proj_nonsq(A, B, param)
% computes the eigenvalues and the corresponding eigenvectors 
% in a circle region for a matrix pencil zB-A
%
% inputs
% A, B: complex matrices
% param.center: radius of the circle
% param.radius: 
% param.L: number of columns of a random matrix, which determines the size of the Rayleigh-Ritz space
% param.M: order of complex moments
% param.N: number of quadrature points for N-point trapezoidal rule
%
% outputs
% X: column vectors are the computed eigenvectors
% lmd: computed eigenvalues
%
% Reference:
% Keiichi Morikuni, 
% Projection method for eigenvalue problems of linear nonsquare matrix pencils,
% SIAM Journal on Matrix Analysis and Applications, Volume 42, Number 3,
% pp. 1381-1400, September 20, 2021.

rng(0, 'twister');

% check the sizes of A, B
[mA, nA] = size(A); 
[mB, nB] = size(B);

if mA ~= mB || nA ~= nB
    error('The sizes of A and B are not the same.');
else
    m = mA; n = nA;
end

A = (A - param.center*B) / param.radius; % shift

% Step 1: Set the values of parameters
L = param.L;
M = param.M;
N = param.N;
LM = L * M;

z = exp(2 * pi * 1i * ((0:N-1)+0.5) / N); % quadrature points
w = z / N; % weights

V = randn(m, L) + 1i*randn(m, L); % pseudo-random matrix

%% Step 2: Form the transformation matrix S
S = zeros(n, LM);
for j = 1:N
    if max(m, n) > 1000 % large case
        tmpmat = globalCGLS(z(j)*B - A, V, param);        
    else % small case
        tmpmat = pinv(z(j)*B - A) * V;        
    end
    S = S + kron(w(j)*z(j).^(0:M-1), tmpmat);
end

%% Step 3. SVD for low-rank approximation
[U_S, SGM_S, ~] = svd(S, 0);
sgm_S = diag(SGM_S);
rank_S = nnz(sgm_S > eps*sgm_S(1));

%% Step 4: Projection onto the range of S
T = randn(m, rank_S) + 1i*randn(m, rank_S);
projA = T(:, 1:rank_S)' * A * U_S(:, 1:rank_S);
projB = T(:, 1:rank_S)' * B * U_S(:, 1:rank_S);

% solution of the projected eigenvalue problem
[V2, lmd] = eig(projA, projB, 'vector');

% Step 5: Preprocessing
X = U_S(:, 1:rank_S) * V2;
X = X ./ vecnorm(X); % normalization

% unshift the computed eigenvalues
lmd = lmd*param.radius + param.center;

end