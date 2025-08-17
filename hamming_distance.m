function dist = hamming_distance(a, b)
%HAMMING_DISTANCE Count elements differing beyond a tolerance.
% The threshold 1e-3 accounts for small numerical precision errors.

dist = sum(abs(a - b) > 1e-3);
