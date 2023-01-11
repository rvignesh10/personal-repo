clc
clear 

%% 
s1 = 'a';
s2 = 't';

for i=1:3
    for j=1:3
        
        temp1 = strcat(s1,num2str(i));
        temp1 = strcat(temp1,num2str(j));

        temp2 = strcat(s2,num2str(i));
        temp2 = strcat(temp2,num2str(j));
        
        A(i,j)= str2sym(temp1);
        T(i,j)= str2sym(temp2);
    end
end

Knum = A*T*A;
disp(Knum(2,3));