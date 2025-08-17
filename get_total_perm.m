function matrix_final = get_total_perm(n_cols)

% random number generator
rng('shuffle');

% max number of possible combinations for matrix B
n_max_comb = 0;
for i = 1:n_cols-1
    n_max_comb = n_max_comb + i;
end
% [ [1,2], [1,3], ..., [col-1, col] ]
matrix_all_perm_indexes = get_all_perm(n_cols);

% finding the resulting matrix
matrix_final = zeros([n_max_comb, n_cols]);
for i = 1:n_max_comb
    % setting 1 in the selected indexes
    matrix_final(i, matrix_all_perm_indexes(i, 1)) = 1;
    matrix_final(i, matrix_all_perm_indexes(i, 2)) = 1;

    % randomly choosing which one will be -1
    negative_index = randi(2);
    matrix_final(i, matrix_all_perm_indexes(i, negative_index)) = -1;
end

% Normalize to ensure unit energy per row
matrix_final = (1/sqrt(2)) * matrix_final;

end

