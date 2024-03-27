function x = gaussianEliminationWithPartialPivoting(A, b, tol)
    
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

    % Perform gaussian elimination with partial pivoting
    for k = 1:n-1
        % Find the pivot index
        [~, j] = max(abs(A(k:n, k)));
        j = j + k - 1;
        
        % Swap the k-th and m-th rows in A and permutation vector p
        A([k j], :) = A([j k], :);
        b([k j],:) = b([j k],:);
        
        % Perform the gaussian elimination
        if  A(k, k)==0
            error('Matrix is singular.');
        else
            lVector = A(k+1:n, k) / A(k, k); % Calculate lower triangular part
        end
        A(k+1:n, k+1:n) = A(k+1:n, k+1:n) - lVector * A(k, k+1:n); % Calculate upper triangular part
        b(k+1:n,:) = b(k+1:n,:) - lVector * b(k,:);
    end
    
    % Check for rank deficiency or near singularity
    diagU = abs(diag(A)); % Get the absolute values of the diagonal of U
    maxDiagU = max(diagU); % Find the largest value
    % Check for values less than tol
    if any(diagU / maxDiagU < tol)
        warning('Matrix is close to singular or badly scaled. Results may be inaccurate.');
    end

    % Backward substitution to solve Ux = y
    % where U is upper triangular. y is altered b while forming U.
    %
    % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
    x = b;
    for i = n:-1:2
        x(i)=x(i)/A(i,i);
        x(1:i-1)= x(1:i-1) - A(1:i-1,i) * x(i);
    end
    x(1)=x(1)/A(1,1);
end
