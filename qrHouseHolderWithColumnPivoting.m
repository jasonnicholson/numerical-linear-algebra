function x = qrHouseHolderWithColumnPivoting(A, b, tol)
  % Solve the least squares problem using QR factorization with pivoting
  % Input: A - Rectangular matrix with more rows than columns
  %        b - Column vector
  %        tol - singularity tolerance similar the reciprocal condition number
  % Output: x - Solution to the least squares problem

  [m, n] = size(A);
  n1 = n + 1;
  nb = size(b,2);
  nnb = n + nb;

  % We could do input checking but this function is not that general
  % purpose.
  %
  % assert(m >= n, "Matrix A must have more rows than columns.");
  % assert(m == size(b,1), "A must have the same number of rows as b.");

  % Combining A and b signifies most of the data needed in cache for processing
  Ab = [A b];

  % Initialize permutation vector
  p = 1:n;

  % Perform QR factorization with column pivoting via Householder reflection
  for k = 1:n
    % Find the column with the largest norm
    [~, maxIdx] = max(sum(Ab(k:m, k:n).^2));
    maxIdx = maxIdx + k - 1;

    % Swap columns k and maxIdx in Ab
    Ab(:, [k, maxIdx]) = Ab(:, [maxIdx, k]);
    p([k, maxIdx]) = p([maxIdx, k]);

    % Find the Householder reflection vector, v
    sigma = -sign(Ab(k,k))*norm(Ab(k:m,k));
    v = [sigma - Ab(k,k); -Ab(k+1:m,k)]; % unnormalized
    vNorm = norm(v);
    if vNorm == 0
      error("A is singular to working precision.");
    else
      v = v/norm(v); % normalized
    end

    % Apply the Householder reflection vector to the submatrix, Ab(k:m,k:nnb)
    Ab(k:m, k:nnb) = Ab(k:m, k:nnb) - 2*v*(v'*Ab(k:m, k:nnb));
  end

  % Check for rank deficiency or near singularity
  diagAb = abs(diag(Ab(1:n,1:n))); % Get the absolute values of the diagonal of U
  % Check for values less than tol
  if any(diagAb / max(diagAb) < tol)
    warning('Matrix is close to singular or badly scaled. Results may be inaccurate.');
  end

  % Perform back substitution to solve R*x = Q'*b for x
  % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
  for i = n:-1:2
    Ab(i,n1:nnb) = Ab(i,n1:nnb) / Ab(i,i);
    Ab(1:i-1,n1:nnb) = Ab(1:i-1,n1:nnb) - Ab(1:i-1,i) * Ab(i,n1:nnb);
  end
  Ab(1,n1:nnb) = Ab(1,n1:nnb)/Ab(1,1);

  % Remove the permutation and output
  x = zeros(n,1);
  x(p) = Ab(1:n,n1:nnb);
end