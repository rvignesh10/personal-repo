% File written by heat1d.C
xa=0; xb=1; kappa=0.1; t=0.2; maxErr= 2.174e-05; cpuTimeStep= 1.591e-02;
Nx=320; dx=  3.125000e-03; numGhost=1; n1a=0; n1b=320; nd1a=-1; nd1b=321;
solutionName='trueNN';
x=[-3.125000000000000e-03 0.000000000000000e+00 3.125000000000000e-03 6.250000000000000e-03 9.375000000000001e-03 1.250000000000000e-02 1.562500000000000e-02 1.875000000000000e-02 ...
2.187500000000000e-02 2.500000000000000e-02 2.812500000000000e-02 3.125000000000000e-02 3.437500000000000e-02 3.750000000000001e-02 4.062500000000000e-02 4.375000000000000e-02 ...
4.687500000000000e-02 5.000000000000000e-02 5.312500000000001e-02 5.625000000000000e-02 5.937500000000000e-02 6.250000000000000e-02 6.562500000000000e-02 6.875000000000001e-02 ...
7.187500000000001e-02 7.500000000000001e-02 7.812500000000000e-02 8.125000000000000e-02 8.437500000000001e-02 8.750000000000001e-02 9.062500000000001e-02 9.375000000000000e-02 ...
9.687500000000000e-02 1.000000000000000e-01 1.031250000000000e-01 1.062500000000000e-01 1.093750000000000e-01 1.125000000000000e-01 1.156250000000000e-01 1.187500000000000e-01 ...
1.218750000000000e-01 1.250000000000000e-01 1.281250000000000e-01 1.312500000000000e-01 1.343750000000000e-01 1.375000000000000e-01 1.406250000000000e-01 1.437500000000000e-01 ...
1.468750000000000e-01 1.500000000000000e-01 1.531250000000000e-01 1.562500000000000e-01 1.593750000000000e-01 1.625000000000000e-01 1.656250000000000e-01 1.687500000000000e-01 ...
1.718750000000000e-01 1.750000000000000e-01 1.781250000000000e-01 1.812500000000000e-01 1.843750000000000e-01 1.875000000000000e-01 1.906250000000000e-01 1.937500000000000e-01 ...
1.968750000000000e-01 2.000000000000000e-01 2.031250000000000e-01 2.062500000000000e-01 2.093750000000000e-01 2.125000000000000e-01 2.156250000000000e-01 2.187500000000000e-01 ...
2.218750000000000e-01 2.250000000000000e-01 2.281250000000000e-01 2.312500000000000e-01 2.343750000000000e-01 2.375000000000000e-01 2.406250000000000e-01 2.437500000000000e-01 ...
2.468750000000000e-01 2.500000000000000e-01 2.531250000000000e-01 2.562500000000000e-01 2.593750000000000e-01 2.625000000000000e-01 2.656250000000000e-01 2.687500000000000e-01 ...
2.718750000000000e-01 2.750000000000000e-01 2.781250000000000e-01 2.812500000000000e-01 2.843750000000000e-01 2.875000000000000e-01 2.906250000000000e-01 2.937500000000000e-01 ...
2.968750000000000e-01 3.000000000000000e-01 3.031250000000000e-01 3.062500000000000e-01 3.093750000000000e-01 3.125000000000000e-01 3.156250000000000e-01 3.187500000000000e-01 ...
3.218750000000000e-01 3.250000000000000e-01 3.281250000000000e-01 3.312500000000000e-01 3.343750000000000e-01 3.375000000000000e-01 3.406250000000000e-01 3.437500000000000e-01 ...
3.468750000000000e-01 3.500000000000000e-01 3.531250000000000e-01 3.562500000000000e-01 3.593750000000000e-01 3.625000000000000e-01 3.656250000000000e-01 3.687500000000000e-01 ...
3.718750000000000e-01 3.750000000000000e-01 3.781250000000000e-01 3.812500000000000e-01 3.843750000000000e-01 3.875000000000000e-01 3.906250000000000e-01 3.937500000000000e-01 ...
3.968750000000000e-01 4.000000000000000e-01 4.031250000000000e-01 4.062500000000000e-01 4.093750000000000e-01 4.125000000000000e-01 4.156250000000000e-01 4.187500000000000e-01 ...
4.218750000000000e-01 4.250000000000000e-01 4.281250000000000e-01 4.312500000000000e-01 4.343750000000000e-01 4.375000000000000e-01 4.406250000000000e-01 4.437500000000000e-01 ...
4.468750000000000e-01 4.500000000000000e-01 4.531250000000000e-01 4.562500000000000e-01 4.593750000000000e-01 4.625000000000000e-01 4.656250000000000e-01 4.687500000000000e-01 ...
4.718750000000000e-01 4.750000000000000e-01 4.781250000000000e-01 4.812500000000000e-01 4.843750000000000e-01 4.875000000000000e-01 4.906250000000000e-01 4.937500000000000e-01 ...
4.968750000000000e-01 5.000000000000000e-01 5.031250000000000e-01 5.062500000000000e-01 5.093750000000000e-01 5.125000000000001e-01 5.156250000000000e-01 5.187500000000000e-01 ...
5.218750000000000e-01 5.250000000000000e-01 5.281250000000001e-01 5.312500000000000e-01 5.343750000000000e-01 5.375000000000000e-01 5.406250000000000e-01 5.437500000000001e-01 ...
5.468750000000000e-01 5.500000000000000e-01 5.531250000000000e-01 5.562500000000000e-01 5.593750000000001e-01 5.625000000000000e-01 5.656250000000000e-01 5.687500000000000e-01 ...
5.718750000000000e-01 5.750000000000001e-01 5.781250000000000e-01 5.812500000000000e-01 5.843750000000000e-01 5.875000000000000e-01 5.906250000000001e-01 5.937500000000000e-01 ...
5.968750000000000e-01 6.000000000000001e-01 6.031250000000000e-01 6.062500000000001e-01 6.093750000000000e-01 6.125000000000000e-01 6.156250000000001e-01 6.187500000000000e-01 ...
6.218750000000001e-01 6.250000000000000e-01 6.281250000000000e-01 6.312500000000001e-01 6.343750000000000e-01 6.375000000000001e-01 6.406250000000000e-01 6.437500000000000e-01 ...
6.468750000000001e-01 6.500000000000000e-01 6.531250000000001e-01 6.562500000000000e-01 6.593750000000000e-01 6.625000000000001e-01 6.656250000000000e-01 6.687500000000001e-01 ...
6.718750000000000e-01 6.750000000000000e-01 6.781250000000001e-01 6.812500000000000e-01 6.843750000000001e-01 6.875000000000000e-01 6.906250000000000e-01 6.937500000000001e-01 ...
6.968750000000000e-01 7.000000000000001e-01 7.031250000000000e-01 7.062500000000000e-01 7.093750000000001e-01 7.125000000000000e-01 7.156250000000001e-01 7.187500000000000e-01 ...
7.218750000000000e-01 7.250000000000001e-01 7.281250000000000e-01 7.312500000000001e-01 7.343750000000000e-01 7.375000000000000e-01 7.406250000000001e-01 7.437500000000000e-01 ...
7.468750000000001e-01 7.500000000000000e-01 7.531250000000000e-01 7.562500000000001e-01 7.593750000000000e-01 7.625000000000001e-01 7.656250000000000e-01 7.687500000000000e-01 ...
7.718750000000001e-01 7.750000000000000e-01 7.781250000000001e-01 7.812500000000000e-01 7.843750000000000e-01 7.875000000000001e-01 7.906250000000000e-01 7.937500000000001e-01 ...
7.968750000000000e-01 8.000000000000000e-01 8.031250000000001e-01 8.062500000000000e-01 8.093750000000001e-01 8.125000000000000e-01 8.156250000000000e-01 8.187500000000001e-01 ...
8.218750000000000e-01 8.250000000000001e-01 8.281250000000000e-01 8.312500000000000e-01 8.343750000000001e-01 8.375000000000000e-01 8.406250000000001e-01 8.437500000000000e-01 ...
8.468750000000000e-01 8.500000000000001e-01 8.531250000000000e-01 8.562500000000001e-01 8.593750000000000e-01 8.625000000000000e-01 8.656250000000001e-01 8.687500000000000e-01 ...
8.718750000000001e-01 8.750000000000000e-01 8.781250000000000e-01 8.812500000000001e-01 8.843750000000000e-01 8.875000000000001e-01 8.906250000000000e-01 8.937500000000000e-01 ...
8.968750000000001e-01 9.000000000000000e-01 9.031250000000001e-01 9.062500000000000e-01 9.093750000000000e-01 9.125000000000001e-01 9.156250000000000e-01 9.187500000000001e-01 ...
9.218750000000000e-01 9.250000000000000e-01 9.281250000000001e-01 9.312500000000000e-01 9.343750000000001e-01 9.375000000000000e-01 9.406250000000000e-01 9.437500000000001e-01 ...
9.468750000000000e-01 9.500000000000001e-01 9.531250000000000e-01 9.562500000000000e-01 9.593750000000001e-01 9.625000000000000e-01 9.656250000000001e-01 9.687500000000000e-01 ...
9.718750000000000e-01 9.750000000000001e-01 9.781250000000000e-01 9.812500000000001e-01 9.843750000000000e-01 9.875000000000000e-01 9.906250000000001e-01 9.937500000000000e-01 ...
9.968750000000001e-01 1.000000000000000e+00 1.003125000000000e+00 ];
u=[1.691728820565208e-01 1.692462827304055e-01 1.691728820565208e-01 1.689527437013632e-01 1.685860586091987e-01 1.680731448364410e-01 1.674144472757751e-01 1.666105372702646e-01 ...
1.656621121177785e-01 1.645699944661671e-01 1.633351315997109e-01 1.619585946174626e-01 1.604415775041940e-01 1.587853960947536e-01 1.569914869327343e-01 1.550614060244394e-01 ...
1.529968274892297e-01 1.507995421074204e-01 1.484714557669889e-01 1.460145878104401e-01 1.434310692832628e-01 1.407231410854972e-01 1.378931520280163e-01 1.349435567952074e-01 ...
1.318769138158202e-01 1.286958830438298e-01 1.254032236512375e-01 1.220017916348119e-01 1.184945373388459e-01 1.148845028960783e-01 1.111748195889996e-01 1.073687051338308e-01 ...
1.034694608895310e-01 9.948046899425465e-02 9.540518943174284e-02 9.124715703019182e-02 8.700997839620338e-02 8.269732878647541e-02 7.831294891994697e-02 7.386064173316206e-02 ...
6.934426908166734e-02 6.476774839030403e-02 6.013504925530058e-02 5.545019000111196e-02 5.071723419499380e-02 4.594028712233288e-02 4.112349222579267e-02 3.627102751136141e-02 ...
3.138710192442088e-02 2.647595169897873e-02 2.154183668323109e-02 1.658903664464273e-02 1.162184755774949e-02 6.644577877903045e-03 1.661544804189997e-03 -3.322929465233232e-03 ...
-8.304521482141548e-03 -1.327891029832669e-02 -1.824178121351422e-02 -2.318882951786118e-02 -2.811576422578818e-02 -3.301831179790737e-02 -3.789221984781786e-02 -4.273326083055325e-02 ...
-4.753723570948241e-02 -5.229997759848234e-02 -5.701735537622440e-02 -6.168527726943865e-02 -6.629969440204847e-02 -7.085660430709677e-02 -7.535205439841783e-02 -7.978214539904320e-02 ...
-8.414303472336845e-02 -8.843093981014639e-02 -9.264214140341651e-02 -9.677298677852417e-02 -1.008198929104321e-01 -1.047793495815753e-01 -1.086479224265640e-01 -1.124222559110938e-01 ...
-1.160990762424790e-01 -1.196751942092852e-01 -1.231475079475964e-01 -1.265130056315192e-01 -1.297687680855899e-01 -1.329119713168167e-01 -1.359398889641641e-01 -1.388498946633517e-01 ...
-1.416394643249178e-01 -1.443061783235719e-01 -1.468477235969366e-01 -1.492618956518582e-01 -1.515466004765470e-01 -1.536998563568873e-01 -1.557197955953425e-01 -1.576046661309644e-01 ...
-1.593528330591007e-01 -1.609627800494839e-01 -1.624331106614701e-01 -1.637625495552881e-01 -1.649499435982472e-01 -1.659942628649450e-01 -1.668946015306073e-01 -1.676501786567848e-01 ...
-1.682603388687262e-01 -1.687245529238388e-01 -1.690424181707448e-01 -1.692136588985339e-01 -1.692381265759107e-01 -1.691157999800278e-01 -1.688467852148946e-01 -1.684313156193438e-01 ...
-1.678697515646380e-01 -1.671625801418896e-01 -1.663104147395656e-01 -1.653139945114459e-01 -1.641741837354932e-01 -1.628919710641931e-01 -1.614684686670140e-01 -1.599049112657298e-01 ...
-1.582026550634435e-01 -1.563631765682400e-01 -1.543880713124872e-01 -1.522790524688993e-01 -1.500379493645594e-01 -1.476667058941923e-01 -1.451673788340636e-01 -1.425421360579665e-01 ...
-1.397932546568453e-01 -1.369231189636851e-01 -1.339342184853824e-01 -1.308291457433886e-01 -1.276105940250013e-01 -1.242813550472528e-01 -1.208443165354216e-01 -1.173024597182696e-01 ...
-1.136588567421740e-01 -1.099166680064002e-01 -1.060791394218254e-01 -1.021495995954902e-01 -9.813145694342211e-02 -9.402819673423277e-02 -8.984337806605555e-02 -8.558063077944451e-02 ...
-8.124365230891202e-02 -7.683620447583733e-02 -7.236211022552637e-02 -6.782525031125419e-02 -6.322955992816548e-02 -5.857902529995321e-02 -5.387768022127597e-02 -4.912960255891319e-02 ...
-4.433891071469264e-02 -3.950976005325890e-02 -3.464633929778060e-02 -2.975286689672322e-02 -2.483358736483865e-02 -1.989276760154519e-02 -1.493469318989162e-02 -9.963664679315289e-03 ...
-4.983993855418704e-03 3.054285377261100e-17 4.983993855418764e-03 9.963664679315351e-03 1.493469318989168e-02 1.989276760154525e-02 2.483358736483871e-02 2.975286689672328e-02 ...
3.464633929778065e-02 3.950976005325896e-02 4.433891071469270e-02 4.912960255891324e-02 5.387768022127604e-02 5.857902529995326e-02 6.322955992816554e-02 6.782525031125426e-02 ...
7.236211022552642e-02 7.683620447583739e-02 8.124365230891208e-02 8.558063077944457e-02 8.984337806605563e-02 9.402819673423284e-02 9.813145694342219e-02 1.021495995954903e-01 ...
1.060791394218254e-01 1.099166680064003e-01 1.136588567421741e-01 1.173024597182697e-01 1.208443165354217e-01 1.242813550472528e-01 1.276105940250014e-01 1.308291457433886e-01 ...
1.339342184853824e-01 1.369231189636851e-01 1.397932546568453e-01 1.425421360579665e-01 1.451673788340637e-01 1.476667058941924e-01 1.500379493645594e-01 1.522790524688993e-01 ...
1.543880713124872e-01 1.563631765682400e-01 1.582026550634435e-01 1.599049112657298e-01 1.614684686670139e-01 1.628919710641931e-01 1.641741837354932e-01 1.653139945114459e-01 ...
1.663104147395656e-01 1.671625801418896e-01 1.678697515646380e-01 1.684313156193438e-01 1.688467852148945e-01 1.691157999800278e-01 1.692381265759107e-01 1.692136588985339e-01 ...
1.690424181707448e-01 1.687245529238388e-01 1.682603388687262e-01 1.676501786567848e-01 1.668946015306073e-01 1.659942628649450e-01 1.649499435982472e-01 1.637625495552881e-01 ...
1.624331106614701e-01 1.609627800494839e-01 1.593528330591007e-01 1.576046661309644e-01 1.557197955953426e-01 1.536998563568873e-01 1.515466004765471e-01 1.492618956518582e-01 ...
1.468477235969366e-01 1.443061783235720e-01 1.416394643249178e-01 1.388498946633517e-01 1.359398889641641e-01 1.329119713168167e-01 1.297687680855899e-01 1.265130056315192e-01 ...
1.231475079475964e-01 1.196751942092852e-01 1.160990762424790e-01 1.124222559110937e-01 1.086479224265640e-01 1.047793495815753e-01 1.008198929104321e-01 9.677298677852411e-02 ...
9.264214140341646e-02 8.843093981014635e-02 8.414303472336844e-02 7.978214539904314e-02 7.535205439841780e-02 7.085660430709673e-02 6.629969440204843e-02 6.168527726943861e-02 ...
5.701735537622435e-02 5.229997759848230e-02 4.753723570948235e-02 4.273326083055320e-02 3.789221984781780e-02 3.301831179790733e-02 2.811576422578812e-02 2.318882951786112e-02 ...
1.824178121351416e-02 1.327891029832663e-02 8.304521482141487e-03 3.322929465233169e-03 -1.661544804190058e-03 -6.644577877903110e-03 -1.162184755774955e-02 -1.658903664464280e-02 ...
-2.154183668323115e-02 -2.647595169897880e-02 -3.138710192442094e-02 -3.627102751136147e-02 -4.112349222579272e-02 -4.594028712233295e-02 -5.071723419499385e-02 -5.545019000111202e-02 ...
-6.013504925530062e-02 -6.476774839030409e-02 -6.934426908166740e-02 -7.386064173316213e-02 -7.831294891994703e-02 -8.269732878647548e-02 -8.700997839620343e-02 -9.124715703019191e-02 ...
-9.540518943174292e-02 -9.948046899425474e-02 -1.034694608895311e-01 -1.073687051338309e-01 -1.111748195889997e-01 -1.148845028960784e-01 -1.184945373388460e-01 -1.220017916348120e-01 ...
-1.254032236512376e-01 -1.286958830438299e-01 -1.318769138158203e-01 -1.349435567952074e-01 -1.378931520280164e-01 -1.407231410854973e-01 -1.434310692832629e-01 -1.460145878104402e-01 ...
-1.484714557669890e-01 -1.507995421074205e-01 -1.529968274892298e-01 -1.550614060244395e-01 -1.569914869327343e-01 -1.587853960947537e-01 -1.604415775041940e-01 -1.619585946174627e-01 ...
-1.633351315997109e-01 -1.645699944661672e-01 -1.656621121177786e-01 -1.666105372702647e-01 -1.674144472757752e-01 -1.680731448364411e-01 -1.685860586091988e-01 -1.689527437013633e-01 ...
-1.691728820565210e-01 -1.692462827304056e-01 -1.691728820565210e-01 ];
err=[2.173081939488148e-05 2.174024795556506e-05 2.173081939488148e-05 2.170254189089847e-05 2.165543997122242e-05 2.158955449101125e-05 2.150494259846118e-05 2.140167768432603e-05 ...
2.127984931884302e-05 2.113956317401267e-05 2.098094093142826e-05 2.080412017747768e-05 2.060925428347250e-05 2.039651227267200e-05 2.016607867401374e-05 1.991815336159365e-05 ...
1.965295138165181e-05 1.937070276585729e-05 1.907165233186743e-05 1.875605947092810e-05 1.842419792294899e-05 1.807635553896339e-05 1.771283403156596e-05 1.733394871312464e-05 ...
1.694002822229236e-05 1.653141423903255e-05 1.610846118815239e-05 1.567153593196921e-05 1.522101745202121e-05 1.475729652033794e-05 1.428077536054889e-05 1.379186729907718e-05 ...
1.329099640642901e-05 1.277859712947298e-05 1.225511391468655e-05 1.172100082245145e-05 1.117672113344071e-05 1.062274694652348e-05 1.005955876955288e-05 9.487645102372836e-06 ...
8.907502013168785e-06 8.319632708119153e-06 7.724547095232621e-06 7.122761341466110e-06 6.514797425640824e-06 5.901182685289365e-06 5.282449359447335e-06 4.659134126870331e-06 ...
4.031777640539579e-06 3.400924058878594e-06 2.767120573482786e-06 2.130916934715017e-06 1.492864974912707e-06 8.535181293717061e-07 2.134309567617401e-07 -4.268413422723333e-07 ...
-1.066743406380336e-06 -1.705720195365011e-06 -2.343217471722980e-06 -2.978682281021062e-06 -3.611563431995203e-06 -4.241311974248239e-06 -4.867381674588979e-06 -5.489229490857389e-06 ...
-6.106316042789401e-06 -6.718106080072124e-06 -7.324068946340060e-06 -7.923679039711936e-06 -8.516416268660530e-06 -9.101766502897632e-06 -9.679222019741428e-06 -1.024828194387541e-05 ...
-1.080845268253900e-05 -1.135924835319421e-05 -1.190019120494048e-05 -1.243081203321680e-05 -1.295065058646565e-05 -1.345925596559237e-05 -1.395618701501380e-05 -1.444101270511142e-05 ...
-1.491331250628591e-05 -1.537267675390089e-05 -1.581870700304974e-05 -1.625101637462113e-05 -1.666922989104978e-05 -1.707298480085138e-05 -1.746193089395636e-05 -1.783573080517395e-05 ...
-1.819406030669870e-05 -1.853660858974991e-05 -1.886307853366940e-05 -1.917318696405195e-05 -1.946666489802436e-05 -1.974325777791362e-05 -2.000272569167692e-05 -2.024484358140018e-05 ...
-2.046940143796071e-05 -2.067620448383921e-05 -2.086507334145563e-05 -2.103584418926652e-05 -2.118836890358531e-05 -2.132251518692911e-05 -2.143816668319246e-05 -2.153522307824142e-05 ...
-2.161360018695952e-05 -2.167323002636646e-05 -2.171406087460641e-05 -2.173605731567499e-05 -2.173920027022182e-05 -2.172348701209007e-05 -2.168893117074955e-05 -2.163556271926887e-05 ...
-2.156342794857703e-05 -2.147258942725114e-05 -2.136312594683067e-05 -2.123513245428222e-05 -2.108871996895587e-05 -2.092401548647717e-05 -2.074116186880403e-05 -2.054031771997208e-05 ...
-2.032165724882230e-05 -2.008537011774508e-05 -1.983166127815851e-05 -1.956075079272366e-05 -1.927287364473007e-05 -1.896827953384459e-05 -1.864723265979953e-05 -1.831001149330915e-05 ...
-1.795690853414213e-05 -1.758823005782199e-05 -1.720429584985290e-05 -1.680543892818235e-05 -1.639200525460568e-05 -1.596435343439174e-05 -1.552285440542184e-05 -1.506789111647948e-05 ...
-1.459985819498231e-05 -1.411916160486192e-05 -1.362621829393666e-05 -1.312145583296816e-05 -1.260531204443418e-05 -1.207823462248142e-05 -1.154068074528474e-05 -1.099311667776825e-05 ...
-1.043601736771656e-05 -9.869866033685190e-06 -9.295153745676201e-06 -8.712378999498843e-06 -8.122047283792070e-06 -7.524670642381419e-06 -6.920767229727028e-06 -6.310860861143608e-06 ...
-5.695480559278626e-06 -5.075160094290339e-06 -4.450437521733297e-06 -3.821854715334427e-06 -3.189956897151368e-06 -2.555292164743555e-06 -1.918411015217742e-06 -1.279865868742955e-06 ...
-6.402105882314813e-07 6.162889791585714e-17 6.402105882041292e-07 1.279865868867639e-06 1.918411015339518e-06 2.555292164715642e-06 3.189956897271409e-06 3.821854715457937e-06 ...
4.450437521842930e-06 5.075160094263547e-06 5.695480559108685e-06 6.310860861257937e-06 6.920767229704933e-06 7.524670642486355e-06 8.122047283903944e-06 8.712378999483687e-06 ...
9.295153745788075e-06 9.869866033797064e-06 1.043601736781904e-05 1.099311667772983e-05 1.154068074516134e-05 1.207823462260717e-05 1.260531204444230e-05 1.312145583309391e-05 ...
1.362621829405750e-05 1.411916160485616e-05 1.459985819510315e-05 1.506789111645984e-05 1.552285440554268e-05 1.596435343440071e-05 1.639200525457622e-05 1.680543892831301e-05 ...
1.720429584990841e-05 1.758823005782114e-05 1.795690853414128e-05 1.831001149332708e-05 1.864723265989262e-05 1.896827953383477e-05 1.927287364465492e-05 1.956075079275141e-05 ...
1.983166127807439e-05 2.008537011780144e-05 2.032165724882230e-05 2.054031771995329e-05 2.074116186872973e-05 2.092401548650493e-05 2.108871996890036e-05 2.123513245428222e-05 ...
2.136312594678413e-05 2.147258942719563e-05 2.156342794857703e-05 2.163556271926887e-05 2.168893117072179e-05 2.172348701209007e-05 2.173920027022182e-05 2.173605731567499e-05 ...
2.171406087460641e-05 2.167323002636646e-05 2.161360018697831e-05 2.153522307824142e-05 2.143816668321125e-05 2.132251518698463e-05 2.118836890360409e-05 2.103584418928530e-05 ...
2.086507334148338e-05 2.067620448384904e-05 2.046940143795089e-05 2.024484358140018e-05 2.000272569180758e-05 1.974325777789483e-05 1.946666489807988e-05 1.917318696406092e-05 ...
1.886307853375351e-05 1.853660858980543e-05 1.819406030672645e-05 1.783573080523031e-05 1.746193089393671e-05 1.707298480084241e-05 1.666922989110614e-05 1.625101637467749e-05 ...
1.581870700301217e-05 1.537267675383556e-05 1.491331250637494e-05 1.444101270521518e-05 1.395618701495743e-05 1.345925596563976e-05 1.295065058645669e-05 1.243081203312372e-05 ...
1.190019120497400e-05 1.135924835317136e-05 1.080845268266603e-05 1.024828194377293e-05 9.679222019723066e-06 9.101766503006299e-06 8.516416268562534e-06 7.923679039717273e-06 ...
7.324068946272699e-06 6.718106079964734e-06 6.106316042811495e-06 5.489229490834543e-06 4.867381674655588e-06 4.241311974101032e-06 3.611563431943374e-06 2.978682281087671e-06 ...
2.343217471561895e-06 1.705720195350144e-06 1.066743406294962e-06 4.268413422976704e-07 -2.134309567721036e-07 -8.535181294602428e-07 -1.492864974886529e-06 -2.130916934917147e-06 ...
-2.767120573527675e-06 -3.400924058817696e-06 -4.031777640693725e-06 -4.659134126918690e-06 -5.282449359697875e-06 -5.901182685452693e-06 -6.514797425680002e-06 -7.122761341406440e-06 ...
-7.724547095346950e-06 -8.319632708179148e-06 -8.907502013111570e-06 -9.487645102226165e-06 -1.005955876959900e-05 -1.062274694648014e-05 -1.117672113353380e-05 -1.172100082251102e-05 ...
-1.225511391463830e-05 -1.277859712959382e-05 -1.329099640647470e-05 -1.379186729925930e-05 -1.428077536066973e-05 -1.475729652035096e-05 -1.522101745219841e-05 -1.567153593210393e-05 ...
-1.610846118821687e-05 -1.653141423919096e-05 -1.694002822239441e-05 -1.733394871315154e-05 -1.771283403153650e-05 -1.807635553909320e-05 -1.842419792300365e-05 -1.875605947094519e-05 ...
-1.907165233197845e-05 -1.937070276593074e-05 -1.965295138180041e-05 -1.991815336173243e-05 -2.016607867407822e-05 -2.039651227284835e-05 -2.060925428352801e-05 -2.080412017764506e-05 ...
-2.098094093153032e-05 -2.113956317404939e-05 -2.127984931896301e-05 -2.140167768440930e-05 -2.150494259857221e-05 -2.158955449110348e-05 -2.165543997129587e-05 -2.170254189100949e-05 ...
-2.173081939497371e-05 -2.174024795567608e-05 -2.173081939499250e-05 ];
