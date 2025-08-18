function [W_SINR, B_seq_SINR, B_obtain,K_SINR, Cz_SINR, Czq_SINR] = sinr_search(B_prime_f, alpha, I_Nr_r, Cn_r, H_r, Cx_r,sigma_n,Nt,Nr)


cx_r=diag(Cx_r);
HCH = H_r*Cx_r*H_r';
Cn_r_const = Cn_r;
HrCx = H_r*Cx_r;
 % Computing the SINR per virtual subchannel per signal component
       
        noise_matix = sigma_n^2 * ones(size(B_prime_f*H_r));
       
        SINR_Matrix_reduced = 1/2*( B_prime_f*H_r ).^2   ./  (  1/2*(B_prime_f*H_r ).^2 *(ones(2*Nt,2*Nt)-eye(2*Nt)) + noise_matix);
        B_seq_SINR=I_Nr_r;  
       
        % Sequential design of the comparator network matrix
        for i_alpha=1:alpha
       
            % MSE per subchannel        
            Cn_r_SINR = B_seq_SINR*Cn_r_const*B_seq_SINR';
            Cz_SINR = B_seq_SINR*HCH*B_seq_SINR'+Cn_r_SINR;
            K_SINR = diag(diag(Cz_SINR).^(-1/2));
            Czqx_SINR = sqrt(2/pi)*K_SINR*B_seq_SINR*HrCx;
            Czq_SINR = 2/pi*(asin(K_SINR*real(Cz_SINR)*K_SINR));
            W_SINR = ((Czq_SINR)^(-1))*Czqx_SINR;
            % MSE per user signal
            MSE_k =  real( diag(W_SINR'*Czq_SINR*W_SINR)  - diag( 2*real(W_SINR'*Czqx_SINR)) +cx_r );    
            % Compute maximum MSE
            [max_MSE , max_MSE_user_index  ] = max(MSE_k);
            % Find the subchannel with max SINR
            [max_SINR , max_SINR_virtual_channel]   = max( SINR_Matrix_reduced(:,max_MSE_user_index));
            %% construction of the comparator matrix
        B_seq_SINR=[B_seq_SINR;B_prime_f(max_SINR_virtual_channel,:) ];
      % Actual the SINR matrix
      SINR_Matrix_reduced(max_SINR_virtual_channel,:)=zeros(1,2*Nt);  
     
     
        end
     
       
        % LRA-MMSE optimized network
       
            Cn_r_SINR = B_seq_SINR*Cn_r_const*B_seq_SINR';
            Cz_SINR = B_seq_SINR*HCH*B_seq_SINR'+Cn_r_SINR;
            K_SINR = diag(diag(Cz_SINR).^(-1/2));
            Czqx_SINR = sqrt(2/pi)*K_SINR*B_seq_SINR*HrCx;
            Czq_SINR = 2/pi*(asin(K_SINR*real(Cz_SINR)*K_SINR));
            W_SINR = ((Czq_SINR)^(-1))*Czqx_SINR;
       

       

            B_obtain=B_seq_SINR(2*Nr+1:1:end,:);


end