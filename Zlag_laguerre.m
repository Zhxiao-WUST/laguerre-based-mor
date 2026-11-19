function Zlag = Zlag_laguerre(A,B,alpha,N)
% ============================================================
% Laguerre coefficient matrix
%
% Zlag = [Z0,Z1,...,Z_{N-1}]
%
% (alpha I - A) Z0 = sqrt(2 alpha) B
%
% (alpha I - A) Zk =
%     sqrt(2 alpha) B
%   - 2 alpha sum_{j=0}^{k-1} Zj
%
% ============================================================

A = sparse(A);
B = sparse(B);

[n,m] = size(B);

I = speye(n);

M = alpha*I - A;

% condition check
if condest(M) > 1e12
    warning('(alpha I - A) is ill-conditioned.');
end

% LU factorization
[L,U,P,Q] = lu(M);

solveM = @(X) Q*(U\(L\(P*X)));

% preallocate
Zlag = zeros(n,N*m);

% ------------------------------------------------------------
% first coefficient
% ------------------------------------------------------------
Z0 = sqrt(2*alpha)*solveM(B);

Zlag(:,1:m)=Z0;

% cumulative sum
Zsum = Z0;

% ------------------------------------------------------------
% recurrence
% ------------------------------------------------------------
for k = 1:N-1

    RHS = sqrt(2*alpha)*B ...
        - 2*alpha*Zsum;

    Znext = solveM(RHS);

    cols = k*m+(1:m);

    Zlag(:,cols)=Znext;

    Zsum = Zsum + Znext;

end

end
