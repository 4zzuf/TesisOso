function [W_SINR, B_seq_SINR, B_obtain,K_SINR, Cz_SINR, Czq_SINR] = sinr_search(B_prime_f, alpha, I_Nr_r, Cn_r, H_r, Cx_r,sigma_n,Nt,Nr)


cx_r=diag(Cx_r);
% Precompute products used throughout the iterations
BfHr   = B_prime_f*H_r;
HrCx   = H_r*Cx_r;
HrCxHr = HrCx*H_r';

% Computing the SINR per virtual subchannel per signal component
noise_matix = sigma_n^2 * ones(size(BfHr));
SINR_Matrix_reduced = 1/2*(BfHr).^2 ./ ( 1/2*(BfHr).^2*(ones(2*Nt,2*Nt)-eye(2*Nt)) + noise_matix);

% Initialize sequential structures
B_seq_SINR   = I_Nr_r;
B_seq_Cn_r   = B_seq_SINR*Cn_r;
B_seq_HrCx   = B_seq_SINR*HrCx;
B_seq_HrCxHr = B_seq_SINR*HrCxHr;

Cn_r_SINR = B_seq_Cn_r*B_seq_SINR';
Cz_signal = B_seq_HrCxHr*B_seq_SINR';
Cz_SINR   = Cz_signal + Cn_r_SINR;

% Sequential design of the comparator network matrix
% parfor cannot be employed here due to the sequential dependence on
% previously selected rows, but it could be used if iterations become
% independent.
for i_alpha=1:alpha

    % MSE per subchannel
    K_SINR = diag(diag(Cz_SINR).^(-1/2));
    Czqx_SINR = sqrt(2/pi)*K_SINR*B_seq_HrCx;
    Czq_SINR = 2/pi*(asin(K_SINR*real(Cz_SINR)*K_SINR));
    W_SINR = (Czq_SINR\Czqx_SINR);
    % MSE per user signal
    MSE_k =  real( diag(W_SINR'*Czq_SINR*W_SINR)  - diag( 2*real(W_SINR'*Czqx_SINR)) +cx_r );
    % Compute maximum MSE
    [max_MSE , max_MSE_user_index  ] = max(MSE_k);
    % Find the subchannel with max SINR
    [max_SINR , max_SINR_virtual_channel]   = max( SINR_Matrix_reduced(:,max_MSE_user_index));
    %% construction of the comparator matrix
    b_new = B_prime_f(max_SINR_virtual_channel,:);

    % Store old values for incremental update
    B_seq_old      = B_seq_SINR;
    B_seq_Cn_old   = B_seq_Cn_r;
    B_seq_HrCx_old = B_seq_HrCx;
    B_seq_HrCxHr_old = B_seq_HrCxHr;

    % Update B matrix and related products
    B_seq_SINR   = [B_seq_old; b_new];
    b_new_Cn     = b_new*Cn_r;
    b_new_HrCx   = b_new*HrCx;
    b_new_HrCxHr = b_new*HrCxHr;

    Cn_r_SINR = [Cn_r_SINR, B_seq_Cn_old*b_new'; b_new_Cn*B_seq_old', b_new_Cn*b_new'];
    Cz_signal = [Cz_signal, B_seq_HrCxHr_old*b_new'; b_new_HrCxHr*B_seq_old', b_new_HrCxHr*b_new'];
    Cz_SINR   = Cz_signal + Cn_r_SINR;

    B_seq_Cn_r   = [B_seq_Cn_old;   b_new_Cn];
    B_seq_HrCx   = [B_seq_HrCx_old; b_new_HrCx];
    B_seq_HrCxHr = [B_seq_HrCxHr_old; b_new_HrCxHr];

    % Update the SINR matrix
    SINR_Matrix_reduced(max_SINR_virtual_channel,:)=zeros(1,2*Nt);

end


% LRA-MMSE optimized network
K_SINR = diag(diag(Cz_SINR).^(-1/2));
Czqx_SINR = sqrt(2/pi)*K_SINR*B_seq_HrCx;
Czq_SINR = 2/pi*(asin(K_SINR*real(Cz_SINR)*K_SINR));
W_SINR = (Czq_SINR\Czqx_SINR);


B_obtain=B_seq_SINR(2*Nr+1:1:end,:);


end