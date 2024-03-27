clc; clear; close all;

% Notes
% https://nhigham.com/2021/04/20/what-is-an-lu-factorization/
% https://www.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox
% https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/

%% Simple single matrix test
rng default;
n = 1000;
A = rand(n);
% A = rand(n);
b = rand(n,1);

results = table;
results.Algorithm = ["\";"LU with Partial Pivoting"; "LU with Rook Pivoting";"LU with Complete Pivoting"];
results.Time(1) = timeit(@() A\b,1);
results.Time(2) = timeit(@() luWithPartialPivoting(A,b),1);
results.Time(3) = timeit(@() luWithRookPivoting(A,b),1);
results.Time(4) = timeit(@() luWithCompletePivoting(A,b),1);

display(results)
x1 = A\b;
x2 = luWithPartialPivoting(A,b);
x3 = luWithRookPivoting(A,b);
x4 = luWithCompletePivoting(A,b);

norm(x1-x2)
norm(x1-x3)
norm(x1-x4)

%% Complex test formed from SVD and spherical cordinates
% https://en.wikipedia.org/wiki/Spherical_coordinate_system
