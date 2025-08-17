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
    
    cont=1;
    aux = zeros(size(new_bit_set,1)-1,size(new_bit_set,2));
    for k=1:size(new_bit_set,1)
        if index~=k
            aux(cont,:) = new_bit_set(k,:);
            cont = cont + 1;
        end
    end
    new_bit_set = aux;
    dist = [];
end
end
