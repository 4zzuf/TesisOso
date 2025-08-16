function dist=hamming_distance(a,b)
dist=0;
for i=1:length(a)
    if abs(a(i)-b(i))>10^(-3)
        dist=dist+1;
    end
end
