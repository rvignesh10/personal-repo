% File written by heat1Dmpi.C
xa=0; xb=1; kappa=0.1; t=0.5; maxErr= 3.387e-16; cpuTimeStep= 6.600e-05;
nx=16; dx=  6.250000e-02; numGhost=1; n1a=0; n1b=16; nd1a=-1; nd1b=17;
solutionName='polyDD';
x=[-6.250000000000000e-02 0.000000000000000e+00 6.250000000000000e-02 1.250000000000000e-01 1.875000000000000e-01 2.500000000000000e-01 3.125000000000000e-01 3.750000000000000e-01 ...
4.375000000000000e-01 5.000000000000000e-01 5.625000000000000e-01 6.250000000000000e-01 6.875000000000000e-01 7.500000000000000e-01 8.125000000000000e-01 8.750000000000000e-01 ...
9.375000000000000e-01 1.000000000000000e+00 1.062500000000000e+00 ];
u=[1.115185546875000e+00 1.150000000000000e+00 1.187060546875000e+00 1.226367187500000e+00 1.267919921875000e+00 1.311718750000000e+00 1.357763671875000e+00 1.406054687500000e+00 ...
1.456591796875000e+00 1.509375000000000e+00 1.564404296875000e+00 1.621679687500000e+00 1.681201171875000e+00 1.742968750000000e+00 1.806982421875000e+00 1.873242187500000e+00 ...
1.941748046875000e+00 2.012500000000000e+00 2.085498046875000e+00 ];
err=[-9.150666335777657e-17 0.000000000000000e+00 4.727121472036799e-17 -1.717376241217039e-16 9.107298248878237e-18 1.457167719820518e-16 1.604619215278547e-17 6.418476861114186e-17 ...
6.808789643208968e-17 2.498001805406602e-16 1.652324110867909e-16 3.642919299551295e-17 3.074797361168891e-16 3.122502256758253e-16 2.727852665973529e-16 1.890848588814720e-16 ...
6.114900252818245e-17 -1.110223024625157e-16 3.387047586844716e-16 ];
