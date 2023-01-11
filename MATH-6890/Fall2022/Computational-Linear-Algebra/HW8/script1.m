syms d e
D = [d/e 0 0; 0 2 0; 0 0 d/e];
A = [1 e e; e 2 e; e e 3];
inv(D)*(A*D)