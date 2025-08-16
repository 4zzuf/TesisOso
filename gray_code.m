function gray_code_data = gray_code(bit_set,modulation_size)

gray_code_data = zeros(size(bit_set));
gray_code_data(1,:) = bit_set(1,:);

new_bit_set=bit_set(2:end,:);

for i=2:modulation_size
    for j=1:size(new_bit_set,1)
        dist(j) = hamming_distance(gray_code_data(i-1,:),new_bit_set(j,:));
    end
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
    clear new_bit_set dist
    new_bit_set = aux;
end
end