function x = gaussianEliminationWithCompletePivoting(A, b, tol)

  % We could check inputs. However, I don't think this algorithm will be used that genericly. Possible checks:

  % Get the size of the matrix A and b vector
  n = size(A, 1);
  nb = size(b,2);
  nnb = n + nb;
  n1 = n+1;

  % Combining A, b, and column pivot. This is most of the data needed in cache for processing
  Ab = [A b; 1:n 0];

  % Perform gaussian elimination with partial pivoting
  for k = 1:n-1
    % Find element that is largest magnitude in the submatrix.
    [~,i] = max(abs(Ab(k:n,k:n)),[],"all");
    subMatrixSize = n-k+1;
    [row,col] = ind2sub([subMatrixSize subMatrixSize],i);
    row = row + k - 1;
    col = col + k - 1;

    % Swap the k-th and m-th rows in A and permutation vector p
    Ab([k row], :) = Ab([row k], :);
    Ab(:, [k col]) = Ab(:, [col k]);

    % Perform the gaussian elimination
    if  Ab(k, k)==0
      error('Matrix is singular.');
    else
      lVector = Ab(k+1:n, k) / Ab(k, k); % Calculate lower triangular part
    end
    Ab(k+1:n, k+1:nnb) = Ab(k+1:n, k+1:nnb) - lVector * Ab(k, k+1:nnb); % Calculate upper triangular part
  end

  % Check for rank deficiency or near singularity
  diagU = abs(diag(Ab(1:n,1:n))); % Get the absolute values of the diagonal of U
  % Check for values less than tol
  if any(diagU / max(diagU) < tol)
    warning('Matrix is close to singular or badly scaled. Results may be inaccurate.');
  end

  % Backward substitution to solve Ux = y
  % where U is upper triangular. y is altered b while forming U.
  %
  % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
  for i = n:-1:2
    Ab(i,n1:nnb) = Ab(i,n1:nnb) / Ab(i,i);
    Ab(1:i-1,n1:nnb) = Ab(1:i-1,n1:nnb) - Ab(1:i-1,i) * Ab(i,n1:nnb);
  end
  Ab(1,n1:nnb) = Ab(1,n1:nnb)/Ab(1,1);

  % Remove permutation and output
  x = zeros(n,1);
  x(Ab(n1,1:n)) = Ab(1:n,n1:nnb);
end
