% File written by heat1d.C
xa=0; xb=1; kappa=0.1; t=0.2; maxErr= 6.127e-03; cpuTimeStep= 1.400e-05;
Nx=20; dx=  5.000000e-02; numGhost=1; n1a=0; n1b=20; nd1a=-1; nd1b=21;
solutionName='trueNN';
x=[-5.000000000000000e-02 0.000000000000000e+00 5.000000000000000e-02 1.000000000000000e-01 1.500000000000000e-01 2.000000000000000e-01 2.500000000000000e-01 3.000000000000000e-01 ...
3.500000000000000e-01 4.000000000000000e-01 4.500000000000000e-01 5.000000000000000e-01 5.500000000000000e-01 6.000000000000001e-01 6.500000000000000e-01 7.000000000000001e-01 ...
7.500000000000000e-01 8.000000000000000e-01 8.500000000000001e-01 9.000000000000000e-01 9.500000000000001e-01 1.000000000000000e+00 1.050000000000000e+00 ];
u=[1.562395105746024e-01 1.753517020730275e-01 1.562395105746024e-01 1.030691444429090e-01 2.743104970768796e-02 -5.418665593313825e-02 -1.239923776284410e-01 -1.667693788999993e-01 ...
-1.731928316410416e-01 -1.418625069696520e-01 -7.960800685431398e-02 -4.183518894777742e-17 7.960800685431384e-02 1.418625069696519e-01 1.731928316410415e-01 1.667693788999992e-01 ...
1.239923776284409e-01 5.418665593313809e-02 -2.743104970768820e-02 -1.030691444429093e-01 -1.562395105746028e-01 -1.753517020730278e-01 -1.562395105746029e-01 ];
err=[5.459339169947893e-03 6.127159590577524e-03 5.459339169947893e-03 3.601454045783834e-03 9.584989327681034e-04 -1.893396440735973e-03 -4.332556095909663e-03 -5.827275054999189e-03 ...
-6.051724088579236e-03 -4.956976236024738e-03 -2.781672244510315e-03 -1.074914480453128e-17 2.781672244510233e-03 4.956976236024598e-03 6.051724088579097e-03 5.827275054999088e-03 ...
4.332556095909473e-03 1.893396440735741e-03 -9.584989327682126e-04 -3.601454045784135e-03 -5.459339169948216e-03 -6.127159590577858e-03 -5.459339169948347e-03 ];
