%% Measure average execution time for several selection strategies (full, random, CVX-based, greedy, SINR)
clear all; close all; clc;
rng('shuffle');

%% Parámetros generales
Nt = 4; Nr = 4;
alpha = 2 * Nr;
SNR_dB = 0;  % Evaluate at a single SNR (0 dB) to reduce computation time
SNR = 10.^(SNR_dB / 10);
sigma_x = 1;
channel_realizations = 4;

%% Matrices base
I_Nr_r = eye(2 * Nr);
Cx_r = (1/2) * sigma_x^2 * eye(2 * Nt);
seed_total_perm = randi(2^32-1);
B_all = get_total_perm(2 * Nr, seed_total_perm);
full = size(B_all, 1);
M_prime_full = 2 * Nr + full;
B_full = [I_Nr_r; B_all];

%% Inicialización de tiempos
time_proposed = zeros(length(SNR), channel_realizations);
time_greedy   = zeros(length(SNR), channel_realizations);
time_sinr     = zeros(length(SNR), channel_realizations);
time_full     = zeros(length(SNR), channel_realizations);
time_random   = zeros(length(SNR), channel_realizations);

for i_channel = 1:channel_realizations
    for i = 1:length(SNR)
        H = (randn(Nr, Nt) + 1i * randn(Nr, Nt)) / sqrt(2);
        H_r = [real(H), -imag(H); imag(H), real(H)];
        sigma_n = sqrt(sigma_x^2 / SNR(i));
        Cn_r = (sigma_n^2 / 2) * I_Nr_r;

        %% Red Full
        A = tic;
        Cz_r_full = B_full * (H_r * Cx_r * H_r') * B_full' + B_full * Cn_r * B_full';
        k_r_full = diag(1 ./ sqrt(diag(Cz_r_full)));
        H_eff_r_q_full = sqrt(2 / pi) * k_r_full * B_full * H_r;
        lambda = (2 / pi) * ((pi / 2 - 1) + (sigma_n^2 / (2 * (Nt * sigma_x^2 / 2 + sigma_n^2 / 2))));
        Time = toc(A);
        time_full(i, i_channel) = Time;

        %% Método propuesto (CVX)
        A = tic;
        cvx_begin quiet sdp
            variable Deltao(M_prime_full);
            maximize(0.5 * log_det(eye(2 * Nt) + 1 / lambda * (sigma_x^2 / 2) * ...
                (H_eff_r_q_full' * diag(Deltao) * H_eff_r_q_full)));
            subject to
                Deltao(1:2 * Nr) == 1;
                0 <= Deltao <= 1;
                sum(Deltao) == 2 * Nr + alpha;
        cvx_end
        Time = toc(A);
        time_proposed(i, i_channel) = Time;

        %% Búsqueda Greedy
        A = tic;
        [~, ~] = greedy_search(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r, full);
        Time = toc(A);
        time_greedy(i, i_channel) = Time;

        %% Selección basada en SINR
        A = tic;
        [~, ~] = sinr_search(B_all, alpha, I_Nr_r, Cn_r, H_r, Cx_r, sigma_n, Nt, Nr);
        Time = toc(A);
        time_sinr(i, i_channel) = Time;

        %% Red Aleatoria
        A = tic;
        seed_rand_alpha = randi(2^32-1);
        B_rand_alpha = get_random_perm(alpha, 2 * Nr, seed_rand_alpha);
        B_rand = [I_Nr_r; B_rand_alpha];
        Cz_r_rand = B_rand * (H_r * Cx_r * H_r') * B_rand' + B_rand * Cn_r * B_rand';
        k_r_rand = diag(1 ./ sqrt(diag(Cz_r_rand)));
        H_eff_r_rand = sqrt(2 / pi) * k_r_rand * B_rand * H_r; %#ok<NASGU>
        Time = toc(A);
        time_random(i, i_channel) = Time;
    end
end

%% Promedios de tiempo
avg_time_full     = mean(time_full(:));
avg_time_random   = mean(time_random(:));
avg_time_proposed = mean(time_proposed(:));
avg_time_greedy   = mean(time_greedy(:));
avg_time_sinr     = mean(time_sinr(:));

fprintf('\n--- Average Execution Times ---\n');
fprintf('Full Network:   %.4f seconds\n', avg_time_full);
fprintf('Random Network: %.4f seconds\n', avg_time_random);
fprintf('Proposed (CVX): %.4f seconds\n', avg_time_proposed);
fprintf('Greedy Search:  %.4f seconds\n', avg_time_greedy);
fprintf('SINR Selection: %.4f seconds\n', avg_time_sinr);
