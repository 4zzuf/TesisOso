function m_comb = get_all_perm(cols)
%GET_ALL_PERM Generate all unique pair combinations
%   m_comb = GET_ALL_PERM(cols) returns all unique combinations
%   of two indices from 1 to cols.

% Generate all combinations of two indices without explicit loops
m_comb = nchoosek(1:cols, 2);

end
