function matrix_final = get_random_perm(n_rows, n_cols, seed)
%GET_RANDOM_PERM Generate random signed pair combinations.
%   MATRIX_FINAL = GET_RANDOM_PERM(N_ROWS, N_COLS, SEED) returns a matrix
%   with N_ROWS randomly selected signed combinations of two columns out of
%   N_COLS columns. When SEED is provided, the random number generator is
%   initialised with this value for reproducibility.

if nargin > 2
    rng(seed);
end

% Max number of possible combinations for matrix B
n_max_comb = n_cols * (n_cols - 1) / 2;

% Generate all combinations and select random rows
matrix_all_perm_indexes = nchoosek(1:n_cols, 2);
array_random_rows_indexes = randperm(n_max_comb, n_rows);
selected_perm = matrix_all_perm_indexes(array_random_rows_indexes, :);

% Preallocate result matrix
matrix_final = zeros(n_rows, n_cols);

% Set ones for the selected indices using vectorisation
rows = repelem((1:n_rows)', 2);
cols = selected_perm(:);
matrix_final(sub2ind([n_rows, n_cols], rows, cols)) = 1;

% Randomly choose which index will be -1 for each row
negative_choice = randi(2, n_rows, 1);
negative_cols = selected_perm(sub2ind(size(selected_perm), (1:n_rows)', negative_choice));
matrix_final(sub2ind([n_rows, n_cols], (1:n_rows)', negative_cols)) = -1;

end

