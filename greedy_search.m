function [final_W, final_B, K, Cz_r] = greedy_search(B_all, alpha, I_nr_r, Cn_r, H_r, Cx_r, n_max_comb)
% Alpha first rows of B_all
B_alpha = B_all(1:alpha,:);
% Computing the MSE for the first alpha rows
new_B = [I_nr_r ; B_alpha];
Cn_rand = new_B*Cn_r*new_B';
K = diag(diag(new_B*H_r*Cx_r*H_r'*new_B'+Cn_rand).^(-1/2));
%Czqx = sqrt(pi/2)*Cx_r*H_r'*new_B'*K;
Czqx = sqrt(pi/2)*(1/2)*Cx_r*H_r'*new_B'*K;
Czq = (asin(transpose(K)*real(new_B*H_r*Cx_r*H_r'*new_B'+Cn_rand)*K) + 1i*asin(transpose(K)*imag(new_B*H_r*Cx_r*H_r'*new_B'+Cn_rand)*K))^(-1);
W = Czqx*Czq;
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
    % Making the permutations in the row that was not freezed with the rows
    % that were not in the first B_alpha
    for j=1:n_max_comb
    	% Checking if the j index is already in B_alpha
        if opt_index_rows(1,j) == 1
            continue;
        end
        B_perm_row = B_all(j,:);
        B_result = [B_frozen(1:i-1,:); B_perm_row; B_frozen(i:end,:)];
        
        % Computing the MSE value
        new_B = [I_nr_r ; B_result];
        Cw_r = new_B*Cn_r*new_B';
        
        %K = diag(diag(new_B*H_r*Cx_r*H_r'*new_B'+Cw_r))^(-1/2);
        K = diag(diag(new_B*H_r*Cx_r*H_r'*new_B'+Cw_r).^(-1/2));
        %Czqx = sqrt(pi/2)*Cx_r*H_r'*new_B'*K;
        Czqx = sqrt(pi/2)*(1/2)*Cx_r*H_r'*new_B'*K;
        Czq = pinv(asin(transpose(K)*real(new_B*H_r*Cx_r*H_r'*new_B'+Cw_r)*K) + 1i*asin(transpose(K)*imag(new_B*H_r*Cx_r*H_r'*new_B'+Cw_r)*K));
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
Cz_r = final_B*H_r*Cx_r*H_r'*final_B' + final_B*Cn_r*final_B';
K = diag(diag(Cz_r).^(-1/2));

end
