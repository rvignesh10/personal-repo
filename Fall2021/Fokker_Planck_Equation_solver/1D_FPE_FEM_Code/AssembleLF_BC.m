function B_out = AssembleLF_BC(B_in,b,e,BoundaryDOF)

B_out = B_in;

for i=1:length(b)
    r = e + i-1;
    bool1 = chkDof(BoundaryDOF,r);
    if ~bool1
        B_out(r,1) = b(i);
    else
        B_out(r,1) = 0;
    end
end

end