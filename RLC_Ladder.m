function [A, B, C, D] = RLC_Ladder(N, R1, R2, L1, C1)

    n_states = 2 * N;
    
    % Precompute constants
    inv_R1C1 = 1/(R1*C1);
    inv_C1 = 1/C1;
    inv_L1 = 1/L1;
    R1_term = R1/(L1*(R1+R2));
    R2_term = R2/(L1*(R1+R2));
    R2_L1 = R2/L1;
    
    % Total non-zero elements in A: 2 (first row) + 6*(N-1) (middle) + 2 (last row)
    n_nz = 6*N - 2;
    I = zeros(n_nz, 1); % Row indices
    J = zeros(n_nz, 1); % Column indices
    V = zeros(n_nz, 1); % Values
    idx = 1; % Current index in I,J,V arrays
    
    % First row (capacitor 1)
    I(idx) = 1; J(idx) = 1; V(idx) = -inv_R1C1; idx = idx + 1;
    I(idx) = 1; J(idx) = 2; V(idx) = -inv_C1;    idx = idx + 1;
    
    % Middle segments (k = 1 to N-1)
    for k = 1:N-1
        % Inductor row (2k)
        row_ind = 2*k;
        I(idx) = row_ind; J(idx) = 2*k-1; V(idx) = inv_L1;    idx = idx + 1;
        I(idx) = row_ind; J(idx) = 2*k;   V(idx) = -R1_term;  idx = idx + 1;
        I(idx) = row_ind; J(idx) = 2*k+1; V(idx) = -R2_term;  idx = idx + 1;
        
        % Capacitor row (2k+1)
        row_cap = 2*k+1;
        I(idx) = row_cap; J(idx) = 2*k;   V(idx) = R2_term;   idx = idx + 1;
        I(idx) = row_cap; J(idx) = 2*k+1; V(idx) = -R2_term;  idx = idx + 1;
        I(idx) = row_cap; J(idx) = 2*k+2; V(idx) = -inv_C1;   idx = idx + 1;
    end
    
    % Last row (inductor N)
    I(idx) = 2*N; J(idx) = 2*N-1; V(idx) = inv_L1;   idx = idx + 1;
    I(idx) = 2*N; J(idx) = 2*N;   V(idx) = -R2_L1;   idx = idx + 1;
    
    % Construct sparse matrices
    X = R1*speye(n_states);
    A = sparse(I, J, V, n_states, n_states);
    B = sparse(1, 1, inv_R1C1, n_states, 1);  % Input affects first state
    C = sparse(1, 1, -1, 1, n_states)/X;     % Output depends on first state
    D = 1/R1;                                 % Direct feedthrough
 
end