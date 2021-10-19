function [A, B] = myrjr2(A, B, ~)
%RJR    Random Jacobi rotation.
%   B = RJR(A) applies a random Jacobi rotation to a square matrix A,
%   preserving its eigenvalues, its singular values, and any symmetry.
%   B = J*A*J', where J is a single random plane rotation.
%
%   B = RJR(A,2) applies two random plane rotations to a rectangular
%   matrix A, preserving its singular values.
%   B = J1*A*J2, where J1 and J2 are two different random plane rotations.

%   R. Schreiber, 1991, CBM, 9-9-93.
%   Copyright 1984-2013 The MathWorks, Inc.

% modified

[m_A, n_A] = size(A);
[m_B, n_B] = size(B);
if nargin == 1 && (m_A~=n_A || m_B~=n_B)
  error(message('MATLAB:rjr:NonSquareMatrix')); 
end

m = m_A;
n = n_A;

if numel(A)<=1 || numel(B)<=1, return, end

% Random rotation of rows.

theta = (2*rand-1)*pi;
c = cos(theta);
s = sin(theta);
i = 1 + fix(rand(1) * m);
j = i;
while (j == i)
   j = 1 + fix(rand(1) * m);
end
A([i j],:) = [c s; -s c] * A([i j],:);
B([i j],:) = [c s; -s c] * B([i j],:);
% Possibly different random rotation of columns.

if nargin > 1
   theta = (2*rand-1)*pi;
   c = cos(theta);
   s = sin(theta);
   i = 1 + fix(rand(1) * n);
   j = i;
   while (j == i)
      j = 1 + fix(rand(1) * n);
   end
end
A(:,[i j]) = A(:,[i j]) * [c -s; s c];
B(:,[i j]) = B(:,[i j]) * [c -s; s c];