function phi = genRBF(fespace)

n = length(fespace); % number of RBFs
%phi = cell(n,n);

for i=1:n
    POE = fespace(i).IntPts; % points to evaluate the RBF in this domain
    InvolvedDOF = fespace(i).Interp.InterpDOF;
    for j=1:length(InvolvedDOF)
        if(InvolvedDOF(j)< 1 || InvolvedDOF(j) > n) % points out of domain
            continue
        else
            k = InvolvedDOF(j);
            if (k>=1 && k<=n)
                involvement = any(fespace(k).Interp.InterpDOF == i); % check involvement of current elemDOF
                if (involvement)
                    centerLoc = fespace(i).Interp.InterpPts(j);
                    count = 1;
                    for o=1:length(fespace(i).Interp.InterpDOF)
                        if (o~=j)
                            zLoc(count) = fespace(i).Interp.InterpPts(o);
                            count = count+1;
                        end
                    end
                    phi(k).x{i} = POE;
                    [v,d]=LocalLagrange(POE,centerLoc,zLoc);
                    phi(k).val{i} = v;
                    phi(k).der{i} = d;
                else
                    phi(k).x{i} = POE;
                    phi(k).val{i} = zeros(length(POE),1);
                    phi(k).der{i} = zeros(length(POE),1);
                end
            end
        end
    end
end

for i=1:n
    stuntN = length(phi(i).x);
    for j=(stuntN+1):n
        phi(i).x{j} = fespace(j).IntPts;
        phi(i).val{j} = zeros(length(fespace(j).IntPts),1);
        phi(i).der{j} = zeros(length(fespace(j).IntPts),1);
    end
    for k=1:n
        if(isempty(phi(i).x{k}))
            phi(i).x{k} = fespace(k).IntPts;
            phi(i).val{k} = zeros(length(fespace(k).IntPts),1);
            phi(i).der{k} = zeros(length(fespace(k).IntPts),1);
        else
            continue
        end
    end
end

end