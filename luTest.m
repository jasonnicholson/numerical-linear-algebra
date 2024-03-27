clc; clear; close all;

% Notes
% https://nhigham.com/2021/04/20/what-is-an-lu-factorization/
% https://www.mathworks.com/matlabcentral/fileexchange/2360-the-matrix-computation-toolbox
% https://nhigham.com/2020/07/14/what-is-the-growth-factor-for-gaussian-elimination/

%% Simple single matrix test
rng default;
n = 1000;
A = rand(n);
% A = rand(5);
b = rand(n,1);

timeit(@() A\b,1)
timeit(@() luWithPartialPivoting(A,b),1)


% norm(x1-x2)
% norm(x1-x3)

%% Complex test formed from SVD and spherical cordinates
% https://en.wikipedia.org/wiki/Spherical_coordinate_system
