function matrix_final = get_total_perm(n_cols, seed)
%GET_TOTAL_PERM Generate all signed pair combinations.
%   MATRIX_FINAL = GET_TOTAL_PERM(N_COLS, SEED) returns a matrix with all
%   possible signed combinations of two columns chosen from N_COLS columns.
%   When SEED is provided, the random number generator is initialised with
%   this value for reproducibility.

if nargin > 1
    rng(seed);
end

% Max number of possible combinations for matrix B
n_max_comb = n_cols * (n_cols - 1) / 2;

% Generate all combinations of two indices
matrix_all_perm_indexes = nchoosek(1:n_cols, 2);

% Preallocate result matrix
matrix_final = zeros(n_max_comb, n_cols);

% Set ones for the selected indices using vectorisation
rows = repelem((1:n_max_comb)', 2);
cols = matrix_all_perm_indexes(:);
matrix_final(sub2ind([n_max_comb, n_cols], rows, cols)) = 1;

% Randomly choose which index will be -1 for each row
negative_choice = randi(2, n_max_comb, 1);
negative_cols = matrix_all_perm_indexes(sub2ind(size(matrix_all_perm_indexes), ...
    (1:n_max_comb)', negative_choice));
matrix_final(sub2ind([n_max_comb, n_cols], (1:n_max_comb)', negative_cols)) = -1;

end

