% File written by heat1dImp.C
xa=0; xb=1; kappa=0.1; t=1; maxErr= 2.076e-02; cpuTimeStep= 1.860e-04;
Nx=20; dx=  5.000000e-02; numGhost=1; n1a=0; n1b=20; nd1a=-1; nd1b=21;
solutionName='trueNN';
x=[-5.000000000000000e-02 0.000000000000000e+00 5.000000000000000e-02 1.000000000000000e-01 1.500000000000000e-01 2.000000000000000e-01 2.500000000000000e-01 3.000000000000000e-01 ...
3.500000000000000e-01 4.000000000000000e-01 4.500000000000000e-01 5.000000000000000e-01 5.500000000000000e-01 6.000000000000001e-01 6.500000000000000e-01 7.000000000000001e-01 ...
7.500000000000000e-01 8.000000000000000e-01 8.500000000000001e-01 9.000000000000000e-01 9.500000000000001e-01 1.000000000000000e+00 1.050000000000000e+00 ];
u=[1.262723206498067e-04 1.417187385522624e-04 1.262723206498067e-04 8.330018449449467e-05 2.216969505156270e-05 -4.379349863408652e-05 -1.002102810515835e-04 -1.347825297813455e-04 ...
-1.399739457120107e-04 -1.146528679102453e-04 -6.433896093787236e-05 -5.644569829741573e-17 6.433896093775896e-05 1.146528679101305e-04 1.399739457118935e-04 1.347825297812254e-04 ...
1.002102810514602e-04 4.379349863395991e-05 -2.216969505169226e-05 -8.330018449462659e-05 -1.262723206499401e-04 -1.417187385523964e-04 -1.262723206499402e-04 ];
err=[2.621322320444682e-06 2.941978817536963e-06 2.621322320444682e-06 1.729251761487129e-06 4.602268824442952e-07 -9.091214517674401e-07 -2.080293172062927e-06 -2.797988125307943e-06 ...
-2.905758176449533e-06 -2.380110860563474e-06 -1.335630433669312e-06 -5.642020542022495e-17 1.335630433555945e-06 2.380110860448575e-06 2.905758176332358e-06 2.797988125187868e-06 ...
2.080293171939559e-06 9.091214516407782e-07 -4.602268825737557e-07 -1.729251761619076e-06 -2.621322320578093e-06 -2.941978817670943e-06 -2.621322320578201e-06 ];