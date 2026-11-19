% ============================================================
% RK-ICOP
% optimal Laguerre parameter
% ============================================================
function alpha = rk_icop(A,B,C,alpha0,N)

tol = 1e-4;
maxit = 50;

m = size(B,2);

a = zeros(maxit+1,1);
a(1)=alpha0;

for i = 1:maxit

    alpha = a(i);

    % Laguerre matrices
    [S_lag,~] = laguerre_SL(alpha,m,N);

    % Laguerre Sylvester
    Zlag = Zlag_laguerre(A,B,alpha,N);

    % reduced system
    Ar = -S_lag';

    % controllability Gramian
    %Xr = lyap(Ar,L_lag'*L_lag);

    % observability-type equation
    Yr = lyap(Ar,eye(N));

    % update alpha
    num = trace(C*A*Zlag*Yr*Zlag'*A'*C');
    den = trace(C*Zlag*Yr*Zlag'*C');

    den = max(real(den),eps);

    a(i+1)=sqrt(real(num/den));

    % stopping criterion
    if abs(a(i+1)-a(i)) <= tol
        break;
    end

end

alpha = a(i+1);

fprintf('Optimal alpha = %.6f\n',alpha);
fprintf('RK-ICOP iterations = %d\n',i);

end
