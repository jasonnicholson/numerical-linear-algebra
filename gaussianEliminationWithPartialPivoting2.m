function x = gaussianEliminationWithPartialPivoting2(A, b, tol)
    
    arguments
        A (:,:) single
        b (:,:) single
        tol (1,1) single = 1e-6;
    end
    
    % We could check inputs. However, I don't think this algorithm will be
    % used that genericly. Possible checks:
    % assert(size(A,1) == size(A,2), "A must be square);
    % assert(size(A,1) == size(b,1), "b must have the number of rows as A")
    % assert(numel(size(A)) == 2, "A must be a matrix")

    % Get the size of the matrix A
    n = size(A, 1);
    m = size(b,2);
    Ab = [A b];

    % Perform gaussian elimination with partial pivoting
    for k = 1:n-1
        % Find the pivot index
        [~, j] = max(abs(Ab(k:n, k)));
        j = j + k - 1;
        
        % Swap the k-th and m-th rows in A and permutation vector p
        Ab([k j], :) = Ab([j k], :);
        
        % Perform the gaussian elimination
        if  Ab(k, k)==0
            error('Matrix is singular.');
        else
            lVector = Ab(k+1:n, k) / Ab(k, k); % Calculate lower triangular part
        end
        Ab(k+1:n, k+1:n+m) = Ab(k+1:n, k+1:n+m) - lVector * Ab(k, k+1:n+m); % Calculate upper triangular part
    end
    
    % Check for rank deficiency or near singularity
    diagU = abs(diag(Ab)); % Get the absolute values of the diagonal of U
    maxDiagU = max(diagU); % Find the largest value
    % Check for values less than tol
    if any(diagU / maxDiagU < tol)
        warning('Matrix is close to singular or badly scaled. Results may be inaccurate.');
    end

    % Backward substitution to solve Ux = y
    % where U is upper triangular. y is altered b while forming U.
    %
    % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
    x = Ab(:,n+1:n+m);
    for i = n:-1:2
        x(i)=x(i)/Ab(i,i);
        x(1:i-1)= x(1:i-1) - Ab(1:i-1,i) * x(i);
    end
    x(1)=x(1)/Ab(1,1);
end
