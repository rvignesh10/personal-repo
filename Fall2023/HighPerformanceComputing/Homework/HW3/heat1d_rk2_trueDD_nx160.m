% File written by heat1d.C
xa=0; xb=1; kappa=0.1; t=0.2; maxErr= 8.706e-05; cpuTimeStep= 2.611e-03;
Nx=160; dx=  6.250000e-03; numGhost=1; n1a=0; n1b=160; nd1a=-1; nd1b=161;
solutionName='trueDD';
x=[-6.250000000000000e-03 0.000000000000000e+00 6.250000000000000e-03 1.250000000000000e-02 1.875000000000000e-02 2.500000000000000e-02 3.125000000000000e-02 3.750000000000001e-02 ...
4.375000000000000e-02 5.000000000000000e-02 5.625000000000000e-02 6.250000000000000e-02 6.875000000000001e-02 7.500000000000001e-02 8.125000000000000e-02 8.750000000000001e-02 ...
9.375000000000000e-02 1.000000000000000e-01 1.062500000000000e-01 1.125000000000000e-01 1.187500000000000e-01 1.250000000000000e-01 1.312500000000000e-01 1.375000000000000e-01 ...
1.437500000000000e-01 1.500000000000000e-01 1.562500000000000e-01 1.625000000000000e-01 1.687500000000000e-01 1.750000000000000e-01 1.812500000000000e-01 1.875000000000000e-01 ...
1.937500000000000e-01 2.000000000000000e-01 2.062500000000000e-01 2.125000000000000e-01 2.187500000000000e-01 2.250000000000000e-01 2.312500000000000e-01 2.375000000000000e-01 ...
2.437500000000000e-01 2.500000000000000e-01 2.562500000000000e-01 2.625000000000000e-01 2.687500000000000e-01 2.750000000000000e-01 2.812500000000000e-01 2.875000000000000e-01 ...
2.937500000000000e-01 3.000000000000000e-01 3.062500000000000e-01 3.125000000000000e-01 3.187500000000000e-01 3.250000000000000e-01 3.312500000000000e-01 3.375000000000000e-01 ...
3.437500000000000e-01 3.500000000000000e-01 3.562500000000000e-01 3.625000000000000e-01 3.687500000000000e-01 3.750000000000000e-01 3.812500000000000e-01 3.875000000000000e-01 ...
3.937500000000000e-01 4.000000000000000e-01 4.062500000000000e-01 4.125000000000000e-01 4.187500000000000e-01 4.250000000000000e-01 4.312500000000000e-01 4.375000000000000e-01 ...
4.437500000000000e-01 4.500000000000000e-01 4.562500000000000e-01 4.625000000000000e-01 4.687500000000000e-01 4.750000000000000e-01 4.812500000000000e-01 4.875000000000000e-01 ...
4.937500000000000e-01 5.000000000000000e-01 5.062500000000000e-01 5.125000000000001e-01 5.187500000000000e-01 5.250000000000000e-01 5.312500000000000e-01 5.375000000000000e-01 ...
5.437500000000001e-01 5.500000000000000e-01 5.562500000000000e-01 5.625000000000000e-01 5.687500000000000e-01 5.750000000000001e-01 5.812500000000000e-01 5.875000000000000e-01 ...
5.937500000000000e-01 6.000000000000001e-01 6.062500000000001e-01 6.125000000000000e-01 6.187500000000000e-01 6.250000000000000e-01 6.312500000000001e-01 6.375000000000001e-01 ...
6.437500000000000e-01 6.500000000000000e-01 6.562500000000000e-01 6.625000000000001e-01 6.687500000000001e-01 6.750000000000000e-01 6.812500000000000e-01 6.875000000000000e-01 ...
6.937500000000001e-01 7.000000000000001e-01 7.062500000000000e-01 7.125000000000000e-01 7.187500000000000e-01 7.250000000000001e-01 7.312500000000001e-01 7.375000000000000e-01 ...
7.437500000000000e-01 7.500000000000000e-01 7.562500000000001e-01 7.625000000000001e-01 7.687500000000000e-01 7.750000000000000e-01 7.812500000000000e-01 7.875000000000001e-01 ...
7.937500000000001e-01 8.000000000000000e-01 8.062500000000000e-01 8.125000000000000e-01 8.187500000000001e-01 8.250000000000001e-01 8.312500000000000e-01 8.375000000000000e-01 ...
8.437500000000000e-01 8.500000000000001e-01 8.562500000000001e-01 8.625000000000000e-01 8.687500000000000e-01 8.750000000000000e-01 8.812500000000001e-01 8.875000000000001e-01 ...
8.937500000000000e-01 9.000000000000000e-01 9.062500000000000e-01 9.125000000000001e-01 9.187500000000001e-01 9.250000000000000e-01 9.312500000000000e-01 9.375000000000000e-01 ...
9.437500000000001e-01 9.500000000000001e-01 9.562500000000000e-01 9.625000000000000e-01 9.687500000000000e-01 9.750000000000001e-01 9.812500000000001e-01 9.875000000000000e-01 ...
9.937500000000000e-01 1.000000000000000e+00 1.006250000000000e+00 ];
u=[-1.000208531211240e-02 0.000000000000000e+00 9.967510215708667e-03 1.990044533501360e-02 2.976435019490120e-02 3.952500908311563e-02 4.914856442492453e-02 5.860163422758459e-02 ...
6.785142787511798e-02 7.686585987173319e-02 8.561366113933924e-02 9.406448748308677e-02 1.021890248486931e-01 1.099590910064367e-01 1.173477333091014e-01 1.243293221847701e-01 ...
1.308796400401619e-01 1.369759652661263e-01 1.425971510538981e-01 1.477236987487172e-01 1.523378254863644e-01 1.564235258779996e-01 1.599666275293274e-01 1.629548402015109e-01 ...
1.653777984433008e-01 1.672270975465024e-01 1.684963227000556e-01 1.691810712416013e-01 1.692789679293463e-01 1.687896731812545e-01 1.677148842529812e-01 1.660583293504681e-01 ...
1.638257546976179e-01 1.610249046039101e-01 1.576654946010982e-01 1.537591777421704e-01 1.493195041794777e-01 1.443618741622407e-01 1.389034846164790e-01 1.329632694926648e-01 ...
1.265618340880207e-01 1.197213835712846e-01 1.124656459578731e-01 1.048197898026242e-01 9.681033689562750e-02 8.846507026397883e-02 7.981293779858166e-02 7.088395184029297e-02 ...
6.170908507372649e-02 5.232016308983529e-02 4.274975398995003e-02 3.303105541421265e-02 2.319777938627827e-02 1.328403537373435e-02 3.324211969874286e-03 -6.647142392754581e-03 ...
-1.659543928324997e-02 -2.648617025199929e-02 -3.628502653288671e-02 -4.595801805284223e-02 -5.547159133591476e-02 -6.479274589288607e-02 -7.388914869269556e-02 -8.272924631859979e-02 ...
-9.128237442002181e-02 -9.951886408042573e-02 -1.074101447322502e-01 -1.149288432619111e-01 -1.220488789611003e-01 -1.287455539950155e-01 -1.349956390737066e-01 -1.407774540293633e-01 ...
-1.460709430200396e-01 -1.508577440989502e-01 -1.551212529080169e-01 -1.588466802747280e-01 -1.620211035125187e-01 -1.646335112467230e-01 -1.666748416106058e-01 -1.681380136789816e-01 ...
-1.690179520303834e-01 -1.693116043525820e-01 -1.690179520303835e-01 -1.681380136789816e-01 -1.666748416106059e-01 -1.646335112467231e-01 -1.620211035125187e-01 -1.588466802747280e-01 ...
-1.551212529080169e-01 -1.508577440989502e-01 -1.460709430200396e-01 -1.407774540293633e-01 -1.349956390737066e-01 -1.287455539950155e-01 -1.220488789611003e-01 -1.149288432619111e-01 ...
-1.074101447322503e-01 -9.951886408042578e-02 -9.128237442002185e-02 -8.272924631859982e-02 -7.388914869269560e-02 -6.479274589288611e-02 -5.547159133591481e-02 -4.595801805284227e-02 ...
-3.628502653288675e-02 -2.648617025199932e-02 -1.659543928325000e-02 -6.647142392754608e-03 3.324211969874260e-03 1.328403537373433e-02 2.319777938627825e-02 3.303105541421263e-02 ...
4.274975398995002e-02 5.232016308983527e-02 6.170908507372649e-02 7.088395184029297e-02 7.981293779858167e-02 8.846507026397882e-02 9.681033689562751e-02 1.048197898026242e-01 ...
1.124656459578731e-01 1.197213835712846e-01 1.265618340880207e-01 1.329632694926648e-01 1.389034846164790e-01 1.443618741622407e-01 1.493195041794777e-01 1.537591777421704e-01 ...
1.576654946010981e-01 1.610249046039101e-01 1.638257546976179e-01 1.660583293504681e-01 1.677148842529812e-01 1.687896731812545e-01 1.692789679293463e-01 1.691810712416013e-01 ...
1.684963227000556e-01 1.672270975465024e-01 1.653777984433008e-01 1.629548402015108e-01 1.599666275293274e-01 1.564235258779996e-01 1.523378254863644e-01 1.477236987487172e-01 ...
1.425971510538981e-01 1.369759652661263e-01 1.308796400401619e-01 1.243293221847701e-01 1.173477333091014e-01 1.099590910064367e-01 1.021890248486931e-01 9.406448748308679e-02 ...
8.561366113933927e-02 7.686585987173322e-02 6.785142787511803e-02 5.860163422758462e-02 4.914856442492458e-02 3.952500908311568e-02 2.976435019490126e-02 1.990044533501366e-02 ...
9.967510215708726e-03 6.217208828649228e-17 -1.000208531211233e-02 ];
err=[-3.970049866579139e-05 0.000000000000000e+00 5.125402262058937e-06 1.023302563315206e-05 1.530515289336211e-05 2.032418995098511e-05 2.527272687249932e-05 3.013359827381956e-05 ...
3.488994286325248e-05 3.952526192957457e-05 4.402347657244036e-05 4.836898347648238e-05 5.254670903587678e-05 5.654216164133166e-05 6.034148194823645e-05 6.393149095165816e-05 ...
6.729973570148578e-05 7.043453249902824e-05 7.332500742490640e-05 7.596113405866232e-05 7.833376825782952e-05 8.043467987733029e-05 8.225658131760524e-05 8.379315280421087e-05 ...
8.503906430912666e-05 8.598999404000624e-05 8.664264343111958e-05 8.699474858555304e-05 8.704508812809655e-05 8.679348744206816e-05 8.624081927466521e-05 8.538900071005233e-05 ...
8.424098651913828e-05 8.280075891006311e-05 8.107331371517819e-05 7.906464306101670e-05 7.678171458352359e-05 7.423244725849894e-05 7.142568393247500e-05 6.837116064898022e-05 ...
6.507947287606067e-05 6.156203875327293e-05 5.783105948448366e-05 5.389947701450629e-05 4.978092913653049e-05 4.548970218561260e-05 4.104068148230681e-05 3.644929969915102e-05 ...
3.173148332756667e-05 2.690359743305310e-05 2.198238888769116e-05 1.698492827945105e-05 1.192855069812556e-05 6.830795603028913e-06 1.709345983369723e-06 -3.418032980901103e-06 ...
-8.533555541984656e-06 -1.361947707936146e-05 -1.865815565220154e-05 -2.363211319464002e-05 -2.852409614425754e-05 -3.331713529019159e-05 -3.799460463582717e-05 -4.254027907089273e-05 ...
-4.693839065236808e-05 -5.117368330057255e-05 -5.523146571850247e-05 -5.909766235306207e-05 -6.275886222034174e-05 -6.620236542462689e-05 -6.941622721237595e-05 -7.238929940500508e-05 ...
-7.511126907040287e-05 -7.757269429571426e-05 -7.976503693906610e-05 -8.168069224710982e-05 -8.331301523363466e-05 -8.465634372993278e-05 -8.570601802556847e-05 -8.645839703165532e-05 ...
-8.691087091116500e-05 -8.706187013204914e-05 -8.691087091119275e-05 -8.645839703165532e-05 -8.570601802557745e-05 -8.465634373000707e-05 -8.331301523364362e-05 -8.168069224711879e-05 ...
-7.976503693910368e-05 -7.757269429569548e-05 -7.511126907045924e-05 -7.238929940499525e-05 -6.941622721231061e-05 -6.620236542473877e-05 -6.275886222043077e-05 -5.909766235310371e-05 ...
-5.523146571841750e-05 -5.117368330068934e-05 -4.693839065248487e-05 -4.254027907099564e-05 -3.799460463587821e-05 -3.331713529011110e-05 -2.852409614439065e-05 -2.363211319468410e-05 ...
-1.865815565225726e-05 -1.361947707939860e-05 -8.533555542026501e-06 -3.418032981016059e-06 1.709345983180918e-06 6.830795602996354e-06 1.192855069816233e-05 1.698492827949823e-05 ...
2.198238888759274e-05 2.690359743288198e-05 3.173148332755728e-05 3.644929969920738e-05 4.104068148237705e-05 4.548970218557994e-05 4.978092913646921e-05 5.389947701448750e-05 ...
5.783105948442239e-05 6.156203875328276e-05 6.507947287596673e-05 6.837116064892385e-05 7.142568393248483e-05 7.423244725847119e-05 7.678171458353341e-05 7.906464306094153e-05 ...
8.107331371507613e-05 8.280075891008189e-05 8.424098651908277e-05 8.538900071009888e-05 8.624081927459092e-05 8.679348744204041e-05 8.704508812809655e-05 8.699474858555304e-05 ...
8.664264343107303e-05 8.598999404002504e-05 8.503906430914546e-05 8.379315280416432e-05 8.225658131767058e-05 8.043467987733029e-05 7.833376825794223e-05 7.596113405874643e-05 ...
7.332500742492520e-05 7.043453249901841e-05 6.729973570141961e-05 6.393149095170961e-05 6.034148194821275e-05 5.654216164125160e-05 5.254670903596090e-05 4.836898347649135e-05 ...
4.402347657261841e-05 3.952526192968686e-05 3.488994286330352e-05 3.013359827376277e-05 2.527272687236002e-05 2.032418995104307e-05 1.530515289333900e-05 1.023302563333204e-05 ...
5.125402262149622e-06 -4.089104698233974e-33 -3.970049866557757e-05 ];