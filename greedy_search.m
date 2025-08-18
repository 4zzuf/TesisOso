function [final_W, final_B, K, Cz_r] = greedy_search(B_all, alpha, I_nr_r, Cn_r, H_r, Cx_r, n_max_comb)
% Precompute constant matrices
HCH = H_r * Cx_r * H_r';
Cn_r_const = Cn_r;
CxHrT_const = sqrt(pi/2)*(1/2) * Cx_r * H_r';

% Alpha first rows of B_all
B_alpha = B_all(1:alpha,:);
% Computing the MSE for the first alpha rows
new_B = [I_nr_r ; B_alpha];
Cn_rand = new_B*Cn_r_const*new_B';
HCH_newB = new_B*HCH*new_B';
K = diag(diag(HCH_newB+Cn_rand).^(-1/2));
Czqx = CxHrT_const*new_B'*K;

% Compute Czq ensuring numerical stability
Czq_matrix = asin(transpose(K) * real(HCH_newB + Cn_rand) * K) ...
    + 1i * asin(transpose(K) * imag(HCH_newB + Cn_rand) * K);

if size(Czq_matrix, 1) == size(Czq_matrix, 2) && cond(Czq_matrix) < 1e10
    Czq = inv(Czq_matrix);
else
    Czq = pinv(Czq_matrix);
end

W = Czqx * Czq;
%equivalent objective MSE
MSE_data = -2*real(trace(Czqx*W')) + trace(W*Czq*W');
%total MSE
MSE_data_2 = trace(Cx_r) - 2*real(trace(Czqx*W')) + trace(W*Czq*W');
% Saving the MSE with the first alpha rows
lower_value = MSE_data;
final_W = W;
final_B = new_B;
lower_value_2 = MSE_data_2;
% Saving the indexes that are already in B_alpha
opt_index_rows = zeros(1,n_max_comb);
opt_index_rows(1,1:alpha) = 1;
for i=1:alpha
    best_index = i;
    % Getting the rows that are not going to be permuted and freeze them
    B_frozen = B_alpha;
    B_frozen(i,:) = [];
    % Precompute base matrices for fixed rows
    row_idx = size(I_nr_r,1) + i;
    Base_rows = [I_nr_r; B_frozen];
    Base_HCH = Base_rows * HCH;
    Base_C = Base_HCH * Base_rows';
    Base_Cn = Base_rows * Cn_r_const * Base_rows';

    % Making the permutations in the row that was not freezed with the rows
    % that were not in the first B_alpha
    for j=1:n_max_comb
        % Checking if the j index is already in B_alpha
        if opt_index_rows(1,j) == 1
            continue;
        end
        B_perm_row = B_all(j,:);
        B_result = [B_frozen(1:i-1,:); B_perm_row; B_frozen(i:end,:)];

        % Update only affected submatrices
        cross_HCH = Base_HCH * B_perm_row';
        self_HCH = B_perm_row * HCH * B_perm_row';
        new_HCH = insert_row_col(Base_C, cross_HCH, self_HCH, row_idx);

        cross_Cn = Base_rows * Cn_r_const * B_perm_row';
        self_Cn = B_perm_row * Cn_r_const * B_perm_row';
        Cw_r = insert_row_col(Base_Cn, cross_Cn, self_Cn, row_idx);

        % Computing the MSE value
        new_B = [I_nr_r ; B_result];
        K = diag(diag(new_HCH + Cw_r).^(-1/2));
        Czqx = CxHrT_const*new_B'*K;
        Czq = pinv(asin(transpose(K)*real(new_HCH + Cw_r)*K) + 1i*asin(transpose(K)*imag(new_HCH + Cw_r)*K));
        W = Czqx*Czq;

        %equivalent objective MSE
        MSE_data = -2*real(trace(Czqx*W')) + trace(W*Czq*W');

        %total MSE
        MSE_data_2 = trace(Cx_r) - 2*real(trace(Czqx*W')) + trace(W*Czq*W');

        % Checking if the curent MSE is the minimum
        if MSE_data < lower_value
            lower_value = MSE_data;
            final_W = W;
            final_B = new_B;
            B_alpha = B_result;

            opt_index_rows(1,best_index) = 0;
            best_index = j;
            opt_index_rows(1,best_index) = 1;
        end

        if MSE_data_2 < lower_value_2
            lower_value_2 = MSE_data_2;
        end

      end

end
% Final Cz_r and K
Cz_r = final_B*HCH*final_B' + final_B*Cn_r_const*final_B';
K = diag(diag(Cz_r).^(-1/2));

end

function updated = insert_row_col(base, cross, self_val, idx)
    n = size(base,1);
    updated = zeros(n+1);
    updated(1:idx-1,1:idx-1) = base(1:idx-1,1:idx-1);
    updated(idx+1:end,idx+1:end) = base(idx:end,idx:end);
    updated(1:idx-1,idx+1:end) = base(1:idx-1,idx:end);
    updated(idx+1:end,1:idx-1) = base(idx:end,1:idx-1);
    updated(1:idx-1,idx) = cross(1:idx-1);
    updated(idx,1:idx-1) = cross(1:idx-1)';
    updated(idx+1:end,idx) = cross(idx:end);
    updated(idx,idx+1:end) = cross(idx:end)';
    updated(idx,idx) = self_val;
end
