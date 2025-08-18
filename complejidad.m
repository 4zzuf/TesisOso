clear all; close all; clc;
%% Parámetros
Nt = 4; Nr = 4;
alpha = 2 * Nr;
SNR_dB = -30:10:40;
SNR = 10.^(SNR_dB / 10);
sigma_x = 1;
channel_realizations = 4;
full = Nr*(2*Nr-1);
M_prime_full = 2 * Nr + full;
%% Identidad y covarianzas
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);
%% Inicialización
time_cvx    = zeros(length(SNR), channel_realizations);
time_greedy = zeros(length(SNR), channel_realizations);
time_sinr   = zeros(length(SNR), channel_realizations);
for i_channel = 1:channel_realizations
    for i = 1:length(SNR)
        H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
        H_r = [real(H), -imag(H); imag(H), real(H)];
        sigma_n = sqrt(sigma_x^2 / SNR(i));
        Cn_r = (sigma_n^2 / 2) * I_Nr_r;
        B_alpha_f = 1/sqrt(2) * randn(full, 2 * Nr);
        B_full = [I_Nr_r; B_alpha_f];
        Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        lambda = (2 / pi)*((pi/2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2/2 + sigma_n^2 / 2))));
        k_r_full = diag(1 ./ sqrt(diag(Cz_r_full)));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;
        %% CVX Optimization
        tic;
        cvx_begin quiet sdp
            variable Deltao(M_prime_full);
            maximize(1/2 * log_det( eye(2*Nt) + 1/lambda *(sigma_x^2/2) * ...
                      (H_eff_r_q_full' * diag(Deltao) * H_eff_r_q_full) ));
            subject to
                Deltao(1:2*Nr) == 1;
                0 <= Deltao <= 1;
                sum(Deltao) == 2*Nr + alpha;
        cvx_end
        time_cvx(i, i_channel) = toc;
        %% Greedy Search
        tic;
        [~, ~] = greedy_search(B_full, alpha, 2*Nr, I_Nr_r, Cn_r, H_r, Cx_r, full, i_channel);
        time_greedy(i, i_channel) = toc;
        %% SINR Selection
        tic;
        [~, ~] = get_B_opt_2_seq_sinr(B_alpha_f, alpha, I_Nr_r, Cn_r, H_r, Cx_r, sigma_n, Nt, Nr);
        time_sinr(i, i_channel) = toc;
    end
end
%% Promedios
avg_time_cvx = mean(time_cvx(:));
avg_time_greedy = mean(time_greedy(:));
avg_time_sinr = mean(time_sinr(:));
fprintf('\n--- Average Execution Times ---\n');
fprintf('CVX Optimization: %.4f seconds\n', avg_time_cvx);
fprintf('Greedy Search:    %.4f seconds\n', avg_time_greedy);
fprintf('SINR Selection:   %.4f seconds\n', avg_time_sinr);