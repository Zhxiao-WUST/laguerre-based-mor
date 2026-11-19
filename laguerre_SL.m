function [S_lag,L_lag] = laguerre_SL(alpha,m,N)
%LAGUERRE_SL Construct the Laguerre generator and input matrices.
%   [S_lag,L_lag] = laguerre_SL(alpha,m,N) returns the matrices defining
%   the finite-dimensional Laguerre realization for scaling parameter
%   alpha, m inputs, and N Laguerre basis functions.
%
%   Inputs:
%       alpha  - positive Laguerre scaling parameter
%       m      - number of system inputs
%       N      - number of Laguerre basis functions
%
%   Outputs:
%       S_lag  - block Laguerre generator matrix
%       L_lag  - Laguerre input matrix
%
%   The construction is consistent with the recursive equations documented
%   below and uses a Kronecker-product structure for the MIMO case.
% ============================================================
% Laguerre generator matrices
%
% Consistent with:
% (alpha I-A)Z0 = sqrt(2alpha)B
%
% (alpha I-A)Zk
% = sqrt(2alpha)B
% -2alpha sum_{j<k}Zj
%
% ============================================================

%% ---------------------------------------------------------
% input matrix
%% ---------------------------------------------------------
% Each Laguerre input block is sqrt(2*alpha) times the m-by-m identity.
L_lag = sqrt(2*alpha) ...
      * repmat(eye(m),1,N);

%% ---------------------------------------------------------
% scalar Laguerre generator
%% ---------------------------------------------------------
% Upper-triangular scalar generator for the N Laguerre basis functions.
S_base = alpha*eye(N) ...
       + 2*alpha*triu(ones(N),1);

%% ---------------------------------------------------------
% MIMO Kronecker structure
%% ---------------------------------------------------------
% Extend the scalar generator to m inputs using a Kronecker product.
S_lag = kron(S_base,eye(m));

end
