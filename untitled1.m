clear all;
close all;
clc;
%% Parámetros generales
Nt = 4; Nr = 12;
alpha = 2 * Nr;
SNR_dB = -30:5:40;
SNR = 10.^(SNR_dB / 10);
sigma_x = 1;
channel_realizations = 4;
full = Nr*(2*Nr-1);
M_prime_full = 2 * Nr + full;
M_prime_random = 2 * Nr + alpha;
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);
% Métricas
I_random = zeros(length(SNR), channel_realizations);
I_select = zeros(length(SNR), channel_realizations);
I_full = zeros(length(SNR), channel_realizations);
I_withoutB = zeros(length(SNR), channel_realizations);
I_greedy = zeros(length(SNR), channel_realizations);
I_sinr = zeros(length(SNR), channel_realizations);
% Modulación QPSK
n_bits = 100;
bits_symbol = 2;
n_symbols = n_bits / bits_symbol;
M = 2^bits_symbol;
constelation_points = exp(1i*pi*(2*(1:M)-1)/M);
constelation_symbols = flipud((dec2bin(M-1:-1:0)-'0'));
gray_code_data = gray_code(constelation_symbols,M);
%% Matrices de comparación
 B_all = (1/sqrt(2)) * get_total_perm(2*Nr);
B_rand = (1/sqrt(2)) * get_random_perm(alpha, 2 * Nr);
B_random_full = [I_Nr_r; B_rand];
B_alpha_f = 1/sqrt(2) * get_random_perm(full, 2 * Nr);
B_full = [I_Nr_r; B_alpha_f];
%% Simulación
h_waitbar = waitbar(0, 'Iniciando simulación...');
for ch = 1:channel_realizations
    waitbar(ch/channel_realizations, h_waitbar, sprintf('Realización del canal %d/%d', ch, channel_realizations));
    H = (randn(Nr,Nt) + 1i*randn(Nr,Nt)) / sqrt(2);
    H_r = [real(H), -imag(H); imag(H), real(H)];
    for i = 1:length(SNR)
        sigma_n = sqrt(sigma_x^2 / SNR(i));
        Cn_r = (sigma_n^2 / 2) * I_Nr_r;
        % ---------- Sin comparadores ----------
        Cz_r = H_r*Cx_r*H_r' + Cn_r;
        K0 = diag(1./sqrt(diag(Cz_r)));
        H_eff_noCN = sqrt(2/pi)*K0*H_r;
        C_eta_noCN = (2/pi)*(asin(K0*Cz_r*K0) - K0*Cz_r*K0) + K0*Cn_r*K0;
        I_withoutB(i,ch) = 0.5*log2(det(eye(2*Nr) + pinv(real(C_eta_noCN)) * ((sigma_x^2/2)*H_eff_noCN*H_eff_noCN')));
        % ---------- Red Random ----------
        Cz_rand = B_random_full * H_r*Cx_r*H_r' * B_random_full' + B_random_full * Cn_r * B_random_full';
        Kr = diag(1 ./ sqrt(diag(Cz_rand)));
        H_eff_rand = sqrt(2/pi)*Kr*B_random_full*H_r;
        C_eta_rand = (2/pi)*(asin(Kr*Cz_rand*Kr) - Kr*Cz_rand*Kr) + Kr*B_random_full*Cn_r*B_random_full'*Kr;
        I_random(i,ch) = 0.5*log2(det(eye(M_prime_random) + pinv(real(C_eta_rand)) * ((sigma_x^2/2)*H_eff_rand*H_eff_rand')));
        % ---------- Red Full ----------
        Cz_full = B_full * H_r*Cx_r*H_r' * B_full' + B_full * Cn_r * B_full';
        Kf = diag(1 ./ sqrt(diag(Cz_full)));
        H_eff_full = sqrt(2/pi)*Kf*B_full*H_r;
        C_eta_full = (2/pi)*(asin(Kf*Cz_full*Kf) - Kf*Cz_full*Kf) + Kf*B_full*Cn_r*B_full'*Kf;
        I_full(i,ch) = 0.5*log2(det(eye(M_prime_full) + pinv(real(C_eta_full)) * ((sigma_x^2/2)*H_eff_full*H_eff_full')));
        % ---------- Red Optimizada (CVX) ----------
        lambda = (2 / pi)*((pi/2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2/2 + sigma_n^2 / 2))));
        cvx_begin quiet sdp
            variable Deltao(M_prime_full)
            maximize(0.5 * log_det( eye(2*Nt) + 1/lambda *(sigma_x^2/2) * ((H_eff_full' * diag(Deltao) * H_eff_full)) ))
            subject to
                for idx = 1:2*Nr
                    Deltao(idx) == 1;
                end
                0 <= Deltao <= 1;
                sum(Deltao) == 2*Nr + alpha;
        cvx_end
        [~, sorted_indices] = maxk(Deltao, 2*Nr + alpha);
        vector_delta_0 = zeros(M_prime_full, 1);
        vector_delta_0(sorted_indices) = 1;
        selected_indices = find(vector_delta_0);
        B_opt = B_full(selected_indices, :);
        Cz_opt = B_opt * H_r * Cx_r * H_r' * B_opt' + B_opt * Cn_r * B_opt';
        Ko = diag(1 ./ sqrt(diag(Cz_opt)));
        H_eff_opt = sqrt(2/pi)*Ko*B_opt*H_r;
        C_eta_opt = (2/pi)*(asin(Ko*Cz_opt*Ko) - Ko*Cz_opt*Ko) + Ko*B_opt*Cn_r*B_opt'*Ko;
        I_select(i,ch) = 0.5*log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_opt)) * ((sigma_x^2/2)*H_eff_opt*H_eff_opt')));
        % ---------- Red Greedy ----------
        [~, B_greedy_norm, K_greedy, Cz_r_greedy] = greedy_search(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r, size(B_all,1));
        H_eff_greedy_norm = sqrt(2/pi)*K_greedy*B_greedy_norm*H_r;
        C_eta_greedy_norm = (2/pi)*(asin(K_greedy*Cz_r_greedy*K_greedy) - K_greedy*Cz_r_greedy*K_greedy) + K_greedy*B_greedy_norm*Cn_r*B_greedy_norm'*K_greedy;
        I_greedy(i,ch) = 0.5*log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_greedy_norm)) * ((sigma_x^2/2)*H_eff_greedy_norm*H_eff_greedy_norm')));
        % ---------- Red SINR ----------
        [~, B_seq_SINR, ~, K_SINR, Cz_SINR, ~] = sinr_search(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r, sigma_n, Nt, Nr);
        H_eff_sinr = sqrt(2/pi) * K_SINR * B_seq_SINR * H_r;
        C_eta_sinr = (2/pi)*(asin(K_SINR*Cz_SINR*K_SINR) - K_SINR*Cz_SINR*K_SINR) + K_SINR*B_seq_SINR*Cn_r*B_seq_SINR'*K_SINR;
        I_sinr(i,ch) = 0.5*log2(det(eye(2*Nr+alpha) + pinv(real(C_eta_sinr)) * ((sigma_x^2/2)*H_eff_sinr*H_eff_sinr')));
    end
end
close(h_waitbar);
% Promedios
I_random_av = mean(I_random, 2);
I_select_av = mean(I_select, 2);
I_full_av = mean(I_full, 2);
I_withoutB_av = mean(I_withoutB, 2);
I_greedy_av = mean(I_greedy, 2);
I_sinr_av = mean(I_sinr, 2);
%% Gráfico final
figure;
plot(SNR_dB, I_random_av, 'g-s', 'LineWidth', 2); hold on;
plot(SNR_dB, I_select_av, 'b-o', 'LineWidth', 2);
plot(SNR_dB, I_full_av, 'r-x', 'LineWidth', 2);
plot(SNR_dB, I_withoutB_av, 'k-.', 'LineWidth', 2);
plot(SNR_dB, I_greedy_av, 'm-d', 'LineWidth', 2);
plot(SNR_dB, I_sinr_av, 'c-^', 'LineWidth', 2);
xlabel('SNR (dB)', 'Interpreter', 'latex');
ylabel('Capacidad (bits/s/Hz)', 'Interpreter', 'latex');
legend({'Random', 'Optimized', 'Full', 'No Comparators', 'Greedy', 'SINR Search'}, ...
    'Interpreter', 'latex', 'Location', 'SouthEast');
title('Capacidad vs. SNR para diferentes topologías de comparadores', 'Interpreter', 'latex');
grid on;



