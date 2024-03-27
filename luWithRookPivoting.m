function x = luWithRookPivoting(A, b, tol)
    
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
    % Initialize the permutation vectors to the identity permutation
    p = 1:n;
    q = 1:n;
    % Initialize LU to the input matrix A
    LU = A;
    
    % Perform LU factorization with partial pivoting
    for k = 1:n-1
        % Find element that is largest magnitude in its row and its column.
        col = k;
        for i=1:n+1-k
            [newMax,row] = max( abs(LU(k:n, col)) );
            row = row + k - 1;
            if i > 1
                if newMax == lastMax
                    break;
                end
            end
            lastMax = newMax;
            [newMax,col] = max( abs(LU(row, k:n)) );
            col = col+k-1;
            if newMax == lastMax
                break;
            end
            lastMax = newMax;
        end
        
        % Permute the Largest element into the pivot position
        LU([k row], :) = LU([row k], :);
        LU(:, [k col]) = LU(:, [col k]);
        p([k row]) = p([row k]);
        q([k col]) = q([col k]);
        
        % Perform the LU factorization
        if  LU(k, k)==0
            error('Matrix is singular.');
        else
            LU(k+1:n, k) = LU(k+1:n, k) / LU(k, k); % Calculate lower triangular part
        end
        LU(k+1:n, k+1:n) = LU(k+1:n, k+1:n) - LU(k+1:n, k) * LU(k, k+1:n); % Calculate upper triangular part 
    end
    
    % Check for rank deficiency or near singularity
    diagU = abs(diag(LU)); % Get the absolute values of the diagonal of U
    maxDiagU = max(diagU); % Find the largest value
    % Check for values less than tol
    if any(diagU / maxDiagU < tol)
        warning('Matrix is close to singular or badly scaled. Results may be inaccurate.');
    end
    
    % Forward substitution to solve Ly = Pb
    % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
    % L(i,i) is assumed to be 1
    % x=y
    x = b(p);
    for i = 1:n-1
        x(i+1:n) = x(i+1:n) - LU(i+1:n,i) * x(i);
    end
    
    % Backward substitution to solve Ux = y
    % Note, this is done a column at a time, hoping to take advantage of matrix stored in column major format.
    for i = n:-1:2
        x(i)=x(i)/LU(i,i);
        x(1:i-1)= x(1:i-1) - LU(1:i-1,i) * x(i);
    end
    x(1)=x(1)/LU(1,1);

    % Undo column permutation
    x(q) = x;
end
