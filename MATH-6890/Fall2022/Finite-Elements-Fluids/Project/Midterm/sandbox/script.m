clc
clear all
%%
choice = 2;
if (choice==1 || choice>=3)
    bc(1)=0; bc(2)=1;
elseif (choice==2)
    bc(1)=0; bc(2)=0;
end
Midterm_Project(bc, choice);