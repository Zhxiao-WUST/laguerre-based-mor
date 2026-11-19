function [H2err, relH2] = relative_H2_error(A,B,C, Ar,Br,Cr)
% ------------------------------------------------------------
% Relative H2 error norm:
%
%   || G - Gr ||_H2 / || G ||_H2
%
% Continuous-time stable LTI systems
%
% INPUT:
%   A,B,C   : full-order system
%   Ar,Br,Cr: reduced-order system
%
% OUTPUT:
%   relH2   : relative H2 error
% ------------------------------------------------------------

%% ---- Check stability ----
if any(real(eig(full(A))) >= 0)
    error('Full system is not asymptotically stable.');
end

if any(real(eig(Ar)) >= 0)
    error('Reduced system is not asymptotically stable.');
end

%% ------------------------------------------------------------
% Full system H2 norm
%% ------------------------------------------------------------
P  = lyap(A, B*B');
H2G_sq = trace(C*P*C');

%% ------------------------------------------------------------
% Reduced system H2 norm
%% ------------------------------------------------------------
Pr = lyap(Ar, Br*Br');
H2Gr_sq = trace(Cr*Pr*Cr');

%% ------------------------------------------------------------
% Cross term  <G,Gr>
%% ------------------------------------------------------------
X = lyap(A, Ar', B*Br');
cross_term = trace(C*X*Cr');

%% ------------------------------------------------------------
% Error H2 norm
%% ------------------------------------------------------------
H2err_sq = H2G_sq + H2Gr_sq - 2*cross_term;
H2err = sqrt(max(H2err_sq,0));

%% ------------------------------------------------------------
% Relative H2 error
%% ------------------------------------------------------------
relH2 = H2err / sqrt(max(H2G_sq,eps));

end
