function gray_code_data = gray_code(bit_set,modulation_size)

gray_code_data = zeros(size(bit_set));
gray_code_data(1,:) = bit_set(1,:);

new_bit_set=bit_set(2:end,:);

for i=2:modulation_size
    % Preassign distance vector to avoid dynamic resizing
    dist = zeros(size(new_bit_set,1),1);

    % Vectorized computation of Hamming distances to all remaining codes
    dist(:) = sum(abs(new_bit_set - gray_code_data(i-1,:)) > 1e-3, 2);
    [~,index]=min(dist);
    gray_code_data(i,:) = new_bit_set(index,:);
    
    new_bit_set(index,:) = [];
    dist = [];
end
end
