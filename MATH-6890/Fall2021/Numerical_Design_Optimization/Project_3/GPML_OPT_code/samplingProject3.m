function samples = samplingProject3(Range,bins,dim)
    
    % use Lattice Hypercube Sampling to generate samples
    t = lhsdesign(bins,dim);
    
    % Scale design variables to their specified range
    for i=1:dim
        t(:,i) = Range(i,1) + (Range(i,2)-Range(i,1)).*t(:,i);
    end

    samples = t;

end