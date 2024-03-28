clear; close all;

%% Simple single matrix test
rng default;
m = 500;
n = 300;
A = rand(m,n);
b = rand(m,1);
tol = 1e-13;

Algorithms = ["\"; "qrHouseholder"; "qrWithColumnPivoting"];

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
