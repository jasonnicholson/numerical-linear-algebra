function x = qrWithColumnPivoting(A, b)
    % Solve the least squares problem using QR factorization with pivoting
    % Input: A - Rectangular matrix with more rows than columns
    %        b - Column vector
    % Output: x - Solution to the least squares problem

    [m, n] = size(A);

    % We could do input checking but this function is not that general
    % purpose.
    %
    % assert(m >= n, "Matrix A must have more rows than columns.");
    % assert(m == size(b,1), "A must have the same number of rows as b.");
    % 
    

    % Initialize permutation vector
    p = 1:n;

    % Perform QR factorization with column pivoting
    for k = 1:n
        % Find the column with the largest norm
        [~, maxIdx] = max(sum(A(k:m, k:n).^2));
        maxIdx = maxIdx + k - 1;

        % Swap columns k and maxIdx in A
        A(:, [k, maxIdx]) = A(:, [maxIdx, k]);
        p([k, maxIdx]) = p([maxIdx, k]);

        % Apply Householder reflection to A(k:m, k)
        sigma = norm(A(k+1:m, k))^2;
        v = [1; A(k+1:m, k)];
        if sigma == 0
            beta = 0;
        else
            alpha = sqrt(A(k, k)^2 + sigma);
            if A(k, k) <= 0
                v(1) = A(k, k) - alpha;
            else
                v(1) = -sigma / (A(k, k) + alpha);
            end
            beta = 2 * v(1)^2 / (sigma + v(1)^2);
            v = v / v(1);
        end
        A(k:m, k:n) = A(k:m, k:n) - beta * v * (v' * A(k:m, k:n));
        A(k+1:m, k) = v(2:end);
    end

    % Apply the transpose of Q to b
    for k = n:-1:1
        v = [1; A(k+1:m, k)];
        beta = 2 / (v' * v);
        b(k:m) = b(k:m) - beta * v * (v' * b(k:m));
    end
    Qb = b(1:n);

    % Perform back substitution to solve Rx = Qb
    x = zeros(n, 1);
    for k = n:-1:1
        x(k) = (Qb(k) - A(k, k+1:n) * x(k+1:n)) / A(k, k);
    end

    % Apply the permutation to x
    x(p) = x;
end