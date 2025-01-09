function phi = genRBF2(fespace)

n = length(fespace);
for i=1:n
    InterpPts = fespace(i).Interp.InterpPts;
    phiCenter = InterpPts(fespace(i).Interp.DOF_ID);
    
    for j=1:length(InterpPts)
        Vtild(:,j) = (InterpPts-phiCenter).^(j-1);
        %Vdtild(:,j) = (j-1)*(InterpPts-phiCenter).^(j-2);
        %Vdtild(:,j) = (InterpPts-phiCenter).^(j-2);
    end

    %Vdtild = Vdtild(:,2:end);
    
    for j=1:length(fespace(i).Interp.InterpPts)
        Vq_k(:,j) = (fespace(i).IntPts - phiCenter).^(j-1);
        %Vdq_k(:,j) = (j-1)*(fespace(i).IntPts - phiCenter).^(j-2);
        %Vdq_k(:,j) = (fespace(i).IntPts - phiCenter).^(j-2);
    end
    
    %Vdq_k = Vdq_k(:,2:end);
    
    Pt     = (Vtild*Vtild')\(Vtild*Vq_k');
    P      = Pt';
    
    
%     Ptd    = (Vdtild*Vdtild')\(Vdtild*Vdq_k');
%     Pd     = Ptd';
    
    for j=1:length(fespace(i).IntPts)
        cLoc = fespace(i).IntPts(j);
        m = 1;
        for k=1:length(fespace(i).IntPts)
            if j==k
                continue
            else
                zLoc(m) = fespace(i).IntPts(k);
                m = m+1;
            end
        end
        POE = fespace(i).IntPts;
        [~,dNdx(:,j)] = LocalLagrange(POE,cLoc,zLoc);
    end
    
    for j=1:length(fespace(i).Interp.InterpDOF)
        phi(fespace(i).Interp.InterpDOF(j)).x{i} = fespace(i).IntPts;
        phi(fespace(i).Interp.InterpDOF(j)).val{i} = P(:,j);
%         phi(fespace(i).Interp.InterpDOF(j)).der{i} = Pd(:,j);
        Pder = P'*dNdx';
        Pder = Pder';
        phi(fespace(i).Interp.InterpDOF(j)).der{i} = Pder(:,j);
        
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