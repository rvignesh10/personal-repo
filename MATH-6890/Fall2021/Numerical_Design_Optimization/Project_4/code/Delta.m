function del = Delta(pt1,pt2,pt3,pt4,L,Nelem)

del = 0;
pt = [pt1 pt2 pt3 pt4];
x = (0:L/Nelem:L)';

for i =1:4
    del = del + pt(i)*cos((2*i-1)*pi*x./(2*L));
end

end