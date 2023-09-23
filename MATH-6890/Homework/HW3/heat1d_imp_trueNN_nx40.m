% File written by heat1dImp.C
xa=0; xb=1; kappa=0.1; t=1; maxErr= 4.742e-03; cpuTimeStep= 3.550e-04;
Nx=40; dx=  2.500000e-02; numGhost=1; n1a=0; n1b=40; nd1a=-1; nd1b=41;
solutionName='trueNN';
x=[-2.500000000000000e-02 0.000000000000000e+00 2.500000000000000e-02 5.000000000000000e-02 7.500000000000001e-02 1.000000000000000e-01 1.250000000000000e-01 1.500000000000000e-01 ...
1.750000000000000e-01 2.000000000000000e-01 2.250000000000000e-01 2.500000000000000e-01 2.750000000000000e-01 3.000000000000000e-01 3.250000000000000e-01 3.500000000000000e-01 ...
3.750000000000000e-01 4.000000000000000e-01 4.250000000000000e-01 4.500000000000000e-01 4.750000000000000e-01 5.000000000000000e-01 5.250000000000000e-01 5.500000000000000e-01 ...
5.750000000000001e-01 6.000000000000001e-01 6.250000000000000e-01 6.500000000000000e-01 6.750000000000000e-01 7.000000000000001e-01 7.250000000000001e-01 7.500000000000000e-01 ...
7.750000000000000e-01 8.000000000000000e-01 8.250000000000001e-01 8.500000000000001e-01 8.750000000000000e-01 9.000000000000000e-01 9.250000000000000e-01 9.500000000000001e-01 ...
9.750000000000001e-01 1.000000000000000e+00 1.025000000000000e+00 ];
u=[1.355852466660022e-04 1.394379277081606e-04 1.355852466660022e-04 1.242401033072815e-04 1.060294320602009e-04 8.195955751709741e-05 5.336058477725894e-05 2.181289762737815e-05 ...
-1.094017371810419e-05 -4.308868932219565e-05 -7.285611709442115e-05 -9.859750423697481e-05 -1.188903775982108e-04 -1.326133497654675e-04 -1.390080871120215e-04 -1.377212154340307e-04 ...
-1.288238474652967e-04 -1.128076531762744e-04 -9.055769001322784e-05 -6.330349448284464e-05 -3.255113776912852e-05 2.215838190018761e-18 3.255113776913291e-05 6.330349448284874e-05 ...
9.055769001323158e-05 1.128076531762776e-04 1.288238474652993e-04 1.377212154340324e-04 1.390080871120224e-04 1.326133497654673e-04 1.188903775982096e-04 9.859750423697262e-05 ...
7.285611709441792e-05 4.308868932219141e-05 1.094017371809897e-05 -2.181289762738426e-05 -5.336058477726585e-05 -8.195955751710504e-05 -1.060294320602090e-04 -1.242401033072901e-04 ...
-1.355852466660110e-04 -1.394379277081695e-04 -1.355852466660111e-04 ];
err=[6.428998496997432e-07 6.611679734351008e-07 6.428998496997432e-07 5.891049779194820e-07 5.027560712734494e-07 3.886247840898720e-07 2.530180294687997e-07 1.034294582597414e-07 ...
-5.187464127720519e-08 -2.043121398765681e-07 -3.454593170970090e-07 -4.675163574542114e-07 -5.637383694660612e-07 -6.288081094299263e-07 -6.591298131561232e-07 -6.530278984695948e-07 ...
-6.108395581474048e-07 -5.348961265925871e-07 -4.293942499235931e-07 -3.001639786415992e-07 -1.543465981119991e-07 2.241331067209537e-18 1.543465981163224e-07 3.001639786457191e-07 ...
4.293942499272659e-07 5.348961265957313e-07 6.108395581499526e-07 6.530278984713295e-07 6.591298131569906e-07 6.288081094297637e-07 5.637383694649363e-07 4.675163574519888e-07 ...
3.454593170937971e-07 2.043121398722720e-07 5.187464127202982e-08 -1.034294582657520e-07 -2.530180294756980e-07 -3.886247840975292e-07 -5.027560712816622e-07 -5.891049779280743e-07 ...
-6.428998497085253e-07 -6.611679734440455e-07 -6.428998497086066e-07 ];
