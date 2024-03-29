clc; clear; close all;

% Notes
% https://nhigham.com/2021/04/20/what-is-an-lu-factorization/
% https://www.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox
% https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/

%% Simple single matrix test
rng default;
n = 800;
A = rand(n);
% A = rand(n);
b = rand(n,1);
tol = 1e-13;

Algorithms = ["\"; "gaussianEliminationWithPartialPivoting"; "gaussianEliminationWithRookPivoting"; ...
  "gaussianEliminationWithCompletePivoting"; "qrHouseholder"; "qrHouseHolderWithColumnNormPivoting"; "qrHouseHolderWithPartialPivoting"];

results = table(Algorithms);
results.Time(1) = timeit(@() A\b,1);
x1 = A\b;
results.NormDiff(1) = 0;
for i=numel(Algorithms):-1:2
    results.Time(i) = timeit(@() feval(Algorithms(i),A,b,tol),1);
    x = feval(Algorithms(i),A,b,tol);
    results.NormDiff(i) = norm(x1 - x);
end

display(results)

%% Complex test formed from SVD and spherical cordinates
% https://en.wikipedia.org/wiki/Spherical_coordinate_system
