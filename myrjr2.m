function [A, B] = myrjr2(A, B, ~)
%myrjr2    Random Jacobi rotation

[m, n] = size(A);

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

% Different random rotation of columns.
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
