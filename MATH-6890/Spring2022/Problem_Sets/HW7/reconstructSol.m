function u = reconstructSol(IDX,v)

for i=1:length(v)
    idx    = IDX{i};
    k      = idx(1);
    j      = idx(2);
    u(k,j) = v(i);
end

end