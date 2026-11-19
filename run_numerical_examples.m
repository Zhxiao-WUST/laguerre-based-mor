   % =========================================================================
% Laguerre-Based Model Order Reduction Framework Demo
% Scenarios: (A) Interpolation, (B) Data-Driven, (C) Min-Phase Preservation
% =========================================================================
clear; clc; close all;

% Add the directory containing this script to the MATLAB search path.
% All MATLAB source files for this example should be placed in this
% same directory; no separate src/ folder is required.
script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

%% 1. Full-order RLC ladder model
% Define the circuit parameters and construct the full-order state-space
% realization. The parameter N determines the size of the ladder model.
 R1 = 0.1;
 R2 = 1;
 L1 = 0.1; 
 C1 = 0.1;
 N = 2000;
 [A, B, C, D] = RLC_Ladder(N, R1, R2, L1, C1);
m = size(B,2);
sys1 = ss(A, B, C, D);
%% 2. Model-based Laguerre MOR
% Select the Laguerre approximation order and an initial scaling parameter.
% The scaling parameter is subsequently updated by rk_icop.
N_A = 20; 
alpha = 1;
alpha = rk_icop(A, B, C, alpha, N_A); 
% Construct the Laguerre generator matrices and the Laguerre coefficient
% matrix Zlag. The Lyapunov-equation solution is used as a consistency
% check for the computed Laguerre coefficient matrix.
[S_lag, L_lag] = laguerre_SL(alpha,m,N_A);
Zlag = Zlag_laguerre(A,B,alpha,N_A);
Zlag_lyap = lyap(A, -S_lag, B*L_lag);
Zlag_err = norm(Zlag-Zlag_lyap,'fro')/norm(Zlag,'fro');
fprintf('Relative error of Zlag: %.3e\n', Zlag_err);
Ar = -S_lag';
Br = L_lag';
Cr = C * Zlag;
sys1_lag = ss(Ar, Br, Cr, D);
%% 3. Projection-based ROM associated with the Laguerre subspace
% Extract an orthonormal basis from Zlag using its thin SVD, and project
% the full-order system matrices onto the resulting subspace.
[V_lag,~,~] = svd(Zlag,0);
Ar_Zlag = V_lag' * A * V_lag;
Br_Zlag = V_lag' * B;
Cr_Zlag = C * V_lag;
sys1_Zlag = ss(Ar_Zlag, Br_Zlag, Cr_Zlag, D);

%% 4. Baseline: balanced-truncation ROM
% Compute a balanced-truncation model of the same target order for
% comparison with the Laguerre-based reduced-order models.
sys1_bt = balred(sys1, N_A);
[Ar_bt, Br_bt, Cr_bt, ~] = ssdata(sys1_bt);
%% ---------------------------------------------------------
% Compare the impulse responses of the full-order model and the ROMs.
% T_end controls the simulation horizon; Nt is the number of time samples.
T_end = 0.06; 
Nt = 2000;
t = linspace(0,T_end,Nt);
[y_full,~] = impulse(sys1,t);
[y_red_lag,~]  = impulse(sys1_lag,t);
[y_red_Zlag,~]  = impulse(sys1_Zlag,t);
[y_red_bt,~]  = impulse(sys1_bt,t);

figure;
plot(t, squeeze(y_full), 'k', 'LineWidth',2); hold on;
plot(t, squeeze(y_red_lag), 'b--', 'LineWidth',1.5);
%plot(t, squeeze(y_red_Zlag), 'm-.', 'LineWidth',1.5);
%plot(t, squeeze(y_red_bt), 'r:', 'LineWidth',1.5);
%legend('Full model','Laguerre ROM','Laguerre Projection','Balanced Truncation','Location','Best');
legend('Full model','Laguerre ROM','Location','Best');
xlabel('Time'); ylabel('Output');
title('Laguerre-based ROM via Z_{lag}');
grid on;

%% ---------------------------------------------------------
% Compute relative errors using the sampled impulse-response vectors.
% The vectorization ensures that the norm is evaluated over all samples.
err_lag = norm(y_full(:) - y_red_lag(:)) / norm(y_full(:));
err_Zlag = norm(y_full(:) - y_red_Zlag(:)) / norm(y_full(:));
err_bt = norm(y_full(:) - y_red_bt(:)) / norm(y_full(:));
fprintf('Relative impulse-response error of Lag: %.3e\n', err_lag);
fprintf('Relative impulse-response error of bt: %.3e\n', err_bt);

%% ---------------------------------------------------------
% Evaluate the residual associated with the Laguerre coefficient
% representation and compare Zlag*Zlag' with the controllability Gramian.
residual = norm(Zlag*L_lag'- B,2)/norm(B,2);
Gram_P = lyap(A, B*B');
res_P = norm(Zlag*Zlag'- Gram_P, 2)/norm(Gram_P,2);
fprintf('Residual norm of Zlag: %.3e\n', residual);
fprintf('Residual norm of Gramian: %.3e\n', res_P);
%% ---------------------------------------------------------
% Evaluate relative H2 errors for the Laguerre-based and balanced-
% truncation reduced-order models using the supplied error routine.
[~, relH2_lag] = relative_H2_error(A,B,C, Ar,Br,Cr);
[~, relH2_bt] = relative_H2_error(A,B,C, Ar_bt,Br_bt,Cr_bt);
%fprintf('H2 error = %.3e\n', H2err);
fprintf('Relative H2 error of lag = %.3e\n', relH2_lag);
fprintf('Relative H2 error of BT = %.3e\n', relH2_bt);
% --- Plot Scenario A ---
figure('Name', 'Scenario A: Interpolation', 'Position', [100, 100, 900, 400], 'Color', 'w');
w = logspace(-1, 3, 1000);
[m_fom] = bode(sys1, w); 
[m_lag] = bode(sys1_lag, w);
 [m_Zlag] = bode(sys1_Zlag, w); 
[m_bt] = bode(sys1_bt, w);

subplot(1, 2, 1);
loglog(w, squeeze(m_fom), 'k', 'LineWidth', 2); hold on;
loglog(w, squeeze(m_lag), 'b--', 'LineWidth', 1.5);
%loglog(w, squeeze(m_Zlag), 'm.-', 'LineWidth', 1.5);
%loglog(w, squeeze(m_bt), 'r:', 'LineWidth', 1.5);
title('Bode Magnitude', 'Interpreter', 'latex');
%legend('FOM', 'Laguerre ROM', 'Laguerre Projection','Balanced Truncation', 'Location', 'best');
legend('FOM', 'Laguerre ROM', 'Location', 'best');
grid on; xlabel('Frequency (rad/s)', 'Interpreter', 'latex');

subplot(1, 2, 2);
err_lag_fre = abs(squeeze(m_fom) - squeeze(m_lag)) ./ squeeze(m_fom);
err_Zlag_fre = abs(squeeze(m_fom) - squeeze(m_Zlag)) ./ squeeze(m_fom);
err_bt_fre  = abs(squeeze(m_fom) - squeeze(m_bt)) ./ squeeze(m_fom);
loglog(w, err_lag_fre, 'b', 'LineWidth', 1.5); hold on;
%loglog(w, err_Zlag_fre, 'm', 'LineWidth', 1.5);
%loglog(w, err_bt_fre, 'r', 'LineWidth', 1.5);
xline(alpha, '--g', 'LineWidth', 1.5);
text(alpha*1.2, min(err_lag_fre)*10, ['$\omega = \alpha = ', num2str(alpha), '$'], 'Interpreter', 'latex', 'Color', 'g');
title('Relative Error $|G(j\omega) - G_r(j\omega)| / |G(j\omega)|$', 'Interpreter', 'latex');
%legend('Laguerre Error','Projection Error', 'BT Error', 'Location', 'best');
legend('Laguerre Error', 'Location', 'best');
grid on; xlabel('Frequency (rad/s)', 'Interpreter', 'latex');

%% =====================================================================
% Scenario B: Data-driven implementation
%
% Reconstruct the Laguerre ROM from sampled impulse-response measurements.
% The coefficient reconstruction below uses the measured time-domain data;
% the full-order matrices are not used in that coefficient calculation.
% =====================================================================

fprintf('\n=========================================\n');
fprintf('Scenario B: Data-driven implementation\n');
fprintf('=========================================\n');

%% ---------------------------------------------------------
% 1. Generate the impulse-response measurement data on the selected
% time grid. The resulting output is converted to a column vector.
t_data = t(:);

[y_meas,~] = impulse(sys1,t_data);

% SISO -> column vector
y_meas = squeeze(y_meas);
y_meas = y_meas(:);

%% ---------------------------------------------------------
% Optional additive measurement noise. Set noise_level to a nonzero
% value to investigate the effect of noisy impulse-response data.
noise_level = 0;

y_meas = y_meas ...
       + noise_level*randn(size(y_meas));

%% ---------------------------------------------------------
% 2. Select Nsamp measurement locations from the available time samples.
% The selected locations are used for plotting and data-driven comparison.
Nsamp = 30;

idx = round(linspace(1,length(t_data),Nsamp));

t_samp = t_data(idx);

y_samp = y_meas(idx);

%% ---------------------------------------------------------
% 3. Compute the Laguerre coefficients from the measured impulse response.
% The coefficient formula is c_k = integral y(t)*phi_k(t) dt, approximated
% numerically here using the trapezoidal rule.
%% ---------------------------------------------------------
c_data = zeros(N_A,1);

for k = 0:N_A-1

    phi = sqrt(2*alpha) ...
        .* exp(-alpha*t_data) ...
        .* laguerreL(k,2*alpha*t_data);

    phi = phi(:);

    c_data(k+1) = trapz( ...
        t_data,...
        y_meas .* phi);

end

%% ---------------------------------------------------------
% 4. Assemble the output matrix of the data-driven ROM from the
% reconstructed Laguerre coefficients. The state and input matrices are
% inherited from the Laguerre realization constructed in Scenario A.
Cr_data = reshape(c_data,1,[]);

sys1_data = ss(Ar,Br,Cr_data,D);

%% ---------------------------------------------------------
% 5. Simulate the impulse response of the data-driven reduced-order model.
[y_data_rom,~] = impulse(sys1_data,t);

%% ---------------------------------------------------------
% 6. Compare the data-driven Laguerre coefficients with the coefficients
% obtained from the model-based Laguerre construction.
coef_err = ...
    norm(Cr-Cr_data,'fro') ...
    / max(norm(Cr,'fro'),eps);

fprintf('Relative coefficient error = %.3e\n',...
    coef_err);

%% ---------------------------------------------------------
% 7. Evaluate the relative H2 error of the data-driven ROM.
[~,relH2_data] = ...
    relative_H2_error( ...
        A,B,C,...
        Ar,Br,Cr_data);

fprintf('Relative H2 error (data-driven) = %.3e\n',...
    relH2_data);

%% ---------------------------------------------------------
% 8. Compare the full-order impulse response, sampled measurements,
% model-based Laguerre ROM, and data-driven Laguerre ROM.
figure( ...
    'Name','Scenario B: Data-driven implementation',...
    'Position',[120 120 920 420],...
    'Color','w');

plot(t,...
    squeeze(y_full),...
    'k',...
    'LineWidth',2);

hold on;

% sampled measurements
plot(...
    t_samp,...
    y_samp,...
    'ro',...
    'MarkerSize',6,...
    'LineWidth',1.2,...
    'MarkerFaceColor','r');

% model-based Laguerre
plot(...
    t,...
    squeeze(y_red_lag),...
    'b--',...
    'LineWidth',1.5);

% data-driven Laguerre
plot(...
    t,...
    squeeze(y_data_rom),...
    'g-.',...
    'LineWidth',1.8);

legend(...
    'Full model',...
    'Sampled measurements',...
    'Model-based Laguerre',...
    'Data-driven Laguerre',...
    'Location','best');

xlabel('Time');

ylabel('Output');

title( ...
    'Scenario B: Data-driven Laguerre reconstruction');

grid on;
box on;

%% ---------------------------------------------------------
% 9. Visualize the model-based and data-driven Laguerre coefficients.
figure( ...
    'Name','Laguerre coefficients',...
    'Position',[160 160 820 380],...
    'Color','w');

stem(0:N_A-1,...
     Cr(:),...
     'b',...
     'LineWidth',1.5);

hold on;

stem(0:N_A-1,...
     Cr_data(:),...
     'r',...
     'LineWidth',1.5);

legend(...
    'Model-based',...
    'Data-driven',...
    'Location','best');

xlabel('Laguerre index');

ylabel('Coefficient');

title('Laguerre coefficients');

grid on;

fprintf('Scenario B completed.\n');

%% =====================================================================
% Scenario C: Minimum-phase preservation
% =====================================================================
% Construct primal and inverse-system Laguerre reduced-order models and
% compare their pole-zero locations and approximation errors.
fprintf('\n==============================\n');
fprintf('Scenario C: Minimum-phase preservation\n');
fprintf('==============================\n');

% ---------------------------------------------------------
% Define the full-order system used in the minimum-phase preservation
% experiment. The matrices are copied from the RLC model constructed above.
% The inverse-dynamics matrix Ax is formed using the direct-feedthrough
% matrix D_mp, as required by the subsequent inverse-system construction.
n_mp = size(A, 1);
A_mp = A;
B_mp = B;
C_mp = C;
D_mp = D;
% Calculate temporary inverse dynamics to check zeros
Ax_temp = A_mp - B_mp*(D_mp\C_mp);

% Recalculate true Ax after shifting
Ax = A_mp - B_mp*(D_mp\C_mp); 
sys_mp =  ss(A_mp, B_mp, C_mp, D_mp);
N_A = 9;
alpha = 14;
% ---------------------------------------------------------
% Standard primal Laguerre construction: compute the Laguerre coefficient
% matrix and form the corresponding reduced-order state-space realization.
m_mp = size(B_mp,2);
[S_mp, L_mp] = laguerre_SL(alpha, m_mp, N_A);

Zp = Zlag_laguerre(A_mp, B_mp, alpha, N_A);
Ar_std = -S_mp';
Br_std = L_mp';
% Note: Assuming Zlag_laguerre projects V. Cr standard is usually -C*V
Cr_std = C_mp * Zp; 
sys_std = ss(Ar_std, Br_std, Cr_std, D_mp);

% ---------------------------------------------------------
% Inverse Laguerre construction:
% First form the state-space matrices of the inverse system, then apply
% the Laguerre construction to that system.
Bx = B_mp / D_mp;
Cx = -(D_mp \ C_mp);

% 2. Project the Inverse System
Zx = Zlag_laguerre(Ax, Bx, alpha, N_A);
Axr_base = -S_mp';
Bxr_base = L_mp';
Cxr_base = Cx * Zx; % Inverse system's reduced C matrix

% Transform the reduced inverse-system realization back to the primal
% input-output representation.
Ar_inv = Axr_base - Bxr_base * D_mp * Cxr_base;
Br_inv = Bxr_base * D_mp;
Cr_inv = -D_mp * Cxr_base;

sys_inv = ss(Ar_inv, Br_inv, Cr_inv, D_mp);

% ---------------------------------------------------------
% Plot the pole-zero maps of the full-order model and both reduced-order
% models. Note that pzmap is applied here to sys_mp as written in this
% script; compatibility depends on the MATLAB model class and release.
figure('Name', 'Scenario C: pole-zero map', 'Color', 'w');
pzmap(sys_mp,'k'); hold on;
pzmap(sys_std,'b');
pzmap(sys_inv,'r');
grid on;
legend('FOM', 'Primal Laguerre', 'Inverse Laguerre', 'Location', 'best');
title('Scenario C: minimum-phase preservation');

% ---------------------------------------------------------
% Count reduced-model transmission zeros in the open right half-plane.
% This provides a numerical check of minimum-phase preservation for the
% two reduced-order models.
z_std = tzero(sys_std);
z_inv = tzero(sys_inv);
nrhp_std = sum(real(z_std)>0);
nrhp_inv = sum(real(z_inv)>0);
fprintf('Primal Laguerre: RHP zeros = %d\n', nrhp_std);
fprintf('Inverse Laguerre: RHP zeros = %d\n', nrhp_inv);
% ---------------------------------------------------------
% Compare the relative H2 and H-infinity approximation errors of the
% primal and inverse Laguerre reduced-order models.
[~, relH2_std] = relative_H2_error(A_mp,B_mp,C_mp, Ar_std,Br_std,Cr_std);
[~, relH2_inv] = relative_H2_error(A_mp,B_mp,C_mp, Ar_inv,Br_inv,Cr_inv);
%fprintf('H2 error = %.3e\n', H2err);
fprintf('Relative H2 error of Primal Laguerre = %.3e\n', relH2_std);
fprintf('Relative H2 error of Inverse Laguerre = %.3e\n', relH2_inv);
% Included for comparison purposes only; no improvements or optimizations 
% have been applied. For large-scale systems, more appropriate algorithms 
% may be selected to enhance computational efficiency.
relH8_std = norm(sys_mp - sys_std, 'inf')/norm(sys_mp,'inf');
relH8_inv = norm(sys_mp - sys_inv, 'inf')/norm(sys_mp,'inf');
fprintf('Relative H_infty error of Primal Laguerre = %.3e\n', relH8_std);
fprintf('Relative H_infty error of Inverse Laguerre = %.3e\n', relH8_inv);
% ---------------------------------------------------------

function y = laguerreL(k,x)
%LAGUERREL Evaluate the (ordinary) Laguerre polynomial of degree k.
%   y = laguerreL(k,x) evaluates L_k(x) elementwise for the entries of x.
%   This local helper uses the finite polynomial expansion and is included
%   so that the main script does not depend on an additional Laguerre
%   polynomial function file.

y = zeros(size(x));

for m = 0:k

    y = y ...
      + (-1)^m ...
      * nchoosek(k,m) ...
      .* x.^m ...
      / factorial(m);

end

end

