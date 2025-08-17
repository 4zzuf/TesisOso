function m_comb = get_all_perm(n_comb, cols)
%GET_ALL_PERM Generate all unique pair combinations
%   m_comb = GET_ALL_PERM(n_comb, cols) returns all unique combinations
%   of two indices from 1 to cols. The input n_comb is kept for
%   compatibility but is not used in the computation.

% Generate all combinations of two indices without explicit loops
m_comb = nchoosek(1:cols, 2);

end
