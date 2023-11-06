% File written by heat1Dmpi.C
xa=0; xb=1; ya=0; yb=1; kappa=0.1; t=0.1; maxErr= 1.349e-15; cpuTimeStep= 1.130e-04;
n1a=0; n1b=16; nd1a=-1; nd1b=17;
n2a=0; n2b=16; nd2a=-1; nd2b=17;
dx(1)=  6.250000e-02; dx(2)=  6.250000e-02; numGhost=1;
option=2; optionName='c-Indexing';
uFinal=[    1.946915473937988e-01,    2.007695312500000e-01,    2.072396430969238e-01,    2.141018829345702e-01,    2.213562507629395e-01,    2.290027465820314e-01,    2.370413703918458e-01,    2.454721221923829e-01,    2.542950019836427e-01,    2.635100097656252e-01,    2.731171455383303e-01,    2.831164093017578e-01,    2.935078010559083e-01,    3.042913208007813e-01,    3.154669685363771e-01,    3.270347442626953e-01,    3.389946479797363e-01,    3.513466796875002e-01,    3.640908393859864e-01;
    1.997636718750000e-01,    2.060000000000000e-01,    2.126386718750000e-01,    2.196796875000000e-01,    2.271230468750000e-01,    2.349687500000000e-01,    2.432167968750000e-01,    2.518671875000000e-01,    2.609199218750000e-01,    2.703750000000000e-01,    2.802324218750000e-01,    2.904921875000000e-01,    3.011542968750001e-01,    3.122187500000000e-01,    3.236855468750000e-01,    3.355546875000000e-01,    3.478261718750000e-01,    3.605000000000000e-01,    3.735761718750000e-01;
    2.071767768859863e-01,    2.136445312500000e-01,    2.205295600891113e-01,    2.278318634033203e-01,    2.355514411926269e-01,    2.436882934570312e-01,    2.522424201965332e-01,    2.612138214111328e-01,    2.706024971008301e-01,    2.804084472656250e-01,    2.906316719055175e-01,    3.012721710205078e-01,    3.123299446105957e-01,    3.238049926757813e-01,    3.356973152160644e-01,    3.480069122314453e-01,    3.607337837219238e-01,    3.738779296875000e-01,    3.874393501281738e-01;
    2.169308624267578e-01,    2.237031250000000e-01,    2.309123077392578e-01,    2.385584106445313e-01,    2.466414337158203e-01,    2.551613769531250e-01,    2.641182403564453e-01,    2.735120239257812e-01,    2.833427276611328e-01,    2.936103515625000e-01,    3.043148956298828e-01,    3.154563598632812e-01,    3.270347442626953e-01,    3.390500488281250e-01,    3.515022735595703e-01,    3.643914184570312e-01,    3.777174835205078e-01,    3.914804687500000e-01,    4.056803741455078e-01;
    2.290259284973146e-01,    2.361757812500000e-01,    2.437869148254394e-01,    2.518593292236328e-01,    2.603930244445801e-01,    2.693880004882813e-01,    2.788442573547363e-01,    2.887617950439453e-01,    2.991406135559082e-01,    3.099807128906250e-01,    3.212820930480957e-01,    3.330447540283203e-01,    3.452686958312988e-01,    3.579539184570312e-01,    3.711004219055176e-01,    3.847082061767578e-01,    3.987772712707519e-01,    4.133076171875000e-01,    4.282992439270023e-01;
    2.434619750976564e-01,    2.510625000000000e-01,    2.591533813476563e-01,    2.677346191406250e-01,    2.768062133789062e-01,    2.863681640625000e-01,    2.964204711914062e-01,    3.069631347656250e-01,    3.179961547851562e-01,    3.295195312500001e-01,    3.415332641601562e-01,    3.540373535156250e-01,    3.670317993164062e-01,    3.805166015625000e-01,    3.944917602539063e-01,    4.089572753906250e-01,    4.239131469726564e-01,    4.393593750000001e-01,    4.552959594726563e-01;
    2.602390022277832e-01,    2.683632812500000e-01,    2.770117073059082e-01,    2.861842803955079e-01,    2.958810005187988e-01,    3.061018676757812e-01,    3.168468818664551e-01,    3.281160430908203e-01,    3.399093513488770e-01,    3.522268066406250e-01,    3.650684089660645e-01,    3.784341583251953e-01,    3.923240547180175e-01,    4.067380981445312e-01,    4.216762886047364e-01,    4.371386260986329e-01,    4.531251106262207e-01,    4.696357421875000e-01,    4.866705207824706e-01;
    2.793570098876954e-01,    2.880781250000000e-01,    2.973618927001953e-01,    3.072083129882812e-01,    3.176173858642578e-01,    3.285891113281250e-01,    3.401234893798828e-01,    3.522205200195312e-01,    3.648802032470703e-01,    3.781025390625000e-01,    3.918875274658203e-01,    4.062351684570312e-01,    4.211454620361328e-01,    4.366184082031250e-01,    4.526540069580078e-01,    4.692522583007813e-01,    4.864131622314453e-01,    5.041367187500001e-01,    5.224229278564456e-01;
    3.008159980773928e-01,    3.102070312500000e-01,    3.202039375305176e-01,    3.308067169189454e-01,    3.420153694152833e-01,    3.538298950195313e-01,    3.662502937316895e-01,    3.792765655517579e-01,    3.929087104797364e-01,    4.071467285156251e-01,    4.219906196594239e-01,    4.374403839111328e-01,    4.534960212707519e-01,    4.701575317382812e-01,    4.874249153137207e-01,    5.052981719970704e-01,    5.237773017883300e-01,    5.428623046875001e-01,    5.625531806945805e-01;
    3.246159667968752e-01,    3.347500000000000e-01,    3.455378417968750e-01,    3.569794921875000e-01,    3.690749511718750e-01,    3.818242187500000e-01,    3.952272949218750e-01,    4.092841796875001e-01,    4.239948730468750e-01,    4.393593750000001e-01,    4.553776855468750e-01,    4.720498046875000e-01,    4.893757324218750e-01,    5.073554687500000e-01,    5.259890136718751e-01,    5.452763671874999e-01,    5.652175292968750e-01,    5.858125000000000e-01,    6.070612792968748e-01;
    3.507569160461425e-01,    3.617070312500000e-01,    3.733636054992676e-01,    3.857266387939454e-01,    3.987961311340332e-01,    4.125720825195313e-01,    4.270544929504395e-01,    4.422433624267578e-01,    4.581386909484864e-01,    4.747404785156250e-01,    4.920487251281739e-01,    5.100634307861328e-01,    5.287845954895020e-01,    5.482122192382813e-01,    5.683463020324707e-01,    5.891868438720703e-01,    6.107338447570800e-01,    6.329873046875001e-01,    6.559472236633305e-01;
    3.792388458251952e-01,    3.910781250000000e-01,    4.036812286376953e-01,    4.170481567382812e-01,    4.311789093017578e-01,    4.460734863281250e-01,    4.617318878173828e-01,    4.781541137695312e-01,    4.953401641845703e-01,    5.132900390625000e-01,    5.320037384033204e-01,    5.514812622070313e-01,    5.717226104736329e-01,    5.927277832031249e-01,    6.144967803955078e-01,    6.370296020507813e-01,    6.603262481689454e-01,    6.843867187500000e-01,    7.092110137939452e-01;
    4.100617561340333e-01,    4.228632812500001e-01,    4.364907112121582e-01,    4.509440460205079e-01,    4.662232856750489e-01,    4.823284301757813e-01,    4.992594795227051e-01,    5.170164337158203e-01,    5.355992927551270e-01,    5.550080566406251e-01,    5.752427253723145e-01,    5.963032989501953e-01,    6.181897773742676e-01,    6.409021606445313e-01,    6.644404487609864e-01,    6.888046417236330e-01,    7.139947395324707e-01,    7.400107421875001e-01,    7.668526496887211e-01;
    4.432256469726561e-01,    4.570625000000000e-01,    4.717920532226563e-01,    4.874143066406250e-01,    5.039292602539063e-01,    5.213369140625000e-01,    5.396372680664063e-01,    5.588303222656251e-01,    5.789160766601563e-01,    5.998945312500000e-01,    6.217656860351563e-01,    6.445295410156250e-01,    6.681860961914063e-01,    6.927353515625000e-01,    7.181773071289064e-01,    7.445119628906250e-01,    7.717393188476562e-01,    7.998593750000000e-01,    8.288721313476562e-01;
    4.787305183410645e-01,    4.936757812500000e-01,    5.095852546691895e-01,    5.264589385986329e-01,    5.442968330383302e-01,    5.630989379882813e-01,    5.828652534484864e-01,    6.035957794189454e-01,    6.252905158996582e-01,    6.479494628906252e-01,    6.715726203918458e-01,    6.961599884033204e-01,    7.217115669250489e-01,    7.482273559570313e-01,    7.757073554992675e-01,    8.041515655517579e-01,    8.335599861145020e-01,    8.639326171875000e-01,    8.952694587707523e-01;
    5.165763702392583e-01,    5.327031250000001e-01,    5.498703155517578e-01,    5.680779418945314e-01,    5.873260040283205e-01,    6.076145019531251e-01,    6.289434356689455e-01,    6.513128051757813e-01,    6.747226104736330e-01,    6.991728515625001e-01,    7.246635284423829e-01,    7.511946411132813e-01,    7.787661895751955e-01,    8.073781738281250e-01,    8.370305938720703e-01,    8.677234497070313e-01,    8.994567413330078e-01,    9.322304687500000e-01,    9.660446319580079e-01;
    5.567632026672361e-01,    5.741445312500000e-01,    5.926472358703614e-01,    6.122713165283203e-01,    6.330167732238771e-01,    6.548836059570314e-01,    6.778718147277834e-01,    7.019813995361329e-01,    7.272123603820801e-01,    7.535646972656251e-01,    7.810384101867676e-01,    8.096334991455079e-01,    8.393499641418457e-01,    8.701878051757814e-01,    9.021470222473146e-01,    9.352276153564454e-01,    9.694295845031738e-01,    1.004752929687500e+00,    1.041197650909424e+00;
    5.992910156249996e-01,    6.180000000000001e-01,    6.379160156250001e-01,    6.590390625000001e-01,    6.813691406250001e-01,    7.049062500000001e-01,    7.296503906250001e-01,    7.556015625000002e-01,    7.827597656250002e-01,    8.111250000000001e-01,    8.406972656250001e-01,    8.714765625000002e-01,    9.034628906250002e-01,    9.366562500000002e-01,    9.710566406250001e-01,    1.006664062500000e+00,    1.043478515625000e+00,    1.081500000000000e+00,    1.120728515625000e+00;
    6.441598091125488e-01,    6.642695312500004e-01,    6.856766548156740e-01,    7.083811798095706e-01,    7.323831062316895e-01,    7.576824340820311e-01,    7.842791633605961e-01,    8.121732940673833e-01,    8.413648262023933e-01,    8.718537597656250e-01,    9.036400947570808e-01,    9.367238311767582e-01,    9.711049690246590e-01,    1.006783508300781e+00,    1.043759449005127e+00,    1.082032791137696e+00,    1.121603534698487e+00,    1.162471679687501e+00,    1.204637226104736e+00;
];
errorG=[   -3.093638762047616e-18,   -6.303985111699718e-18,   -4.144338748590139e-17,   -1.648556643843793e-16,    3.979281842730019e-17,    1.014899427519145e-16,    7.491419255225250e-17,    4.333229467520100e-17,    6.308806762048668e-17,    1.341815113881095e-16,    1.439249889786161e-16,   -4.396137588073140e-17,    1.239685945857074e-16,    3.709223683409646e-17,    8.398760948324056e-17,   -6.841219485440725e-17,   -2.986378302310445e-17,    1.407910246720157e-16,    2.793687454765921e-18;
   -3.514376289981414e-17,   -3.330669073875473e-19,    6.989200884710556e-18,   -1.317695952351983e-17,   -5.320396900820868e-18,    2.803313137178520e-18,    1.202683785894720e-17,    2.151750999601632e-17,   -2.506848895134084e-17,    1.354472090042691e-17,   -3.918740332231608e-18,   -2.028238688112083e-17,    1.829959794807934e-17,    2.470246229790972e-18,   -1.225929080472810e-17,   -2.755434769241560e-17,    1.376156133492401e-17,   -9.992007221626418e-19,   -1.466014809548000e-17;
   -2.589863883789047e-17,    8.073403057196060e-18,   -8.202941499885555e-18,    9.371721606220295e-18,   -2.246933447137313e-17,   -2.129205015424795e-17,   -1.401933378956423e-17,   -6.511853773219497e-19,   -9.775847801618900e-18,   -6.998156394655286e-17,   -6.691536227573271e-17,   -5.941906309429174e-17,    1.134915390290335e-17,    3.103631717946143e-17,   -5.586872449587534e-17,   -2.399069712420016e-17,   -4.319372347316199e-17,    8.751679936302995e-19,   -6.164814549147225e-17;
   -3.170208507600331e-17,   -8.840150833577809e-18,   -2.920420524334111e-18,   -1.311022687980334e-17,   -1.165399428435659e-17,    6.156099935372691e-19,   -3.181256527737958e-17,   -2.317379144481380e-17,    2.320126241735915e-17,   -5.922085738463423e-17,   -4.506487685188715e-17,   -4.868376752079073e-17,   -1.123570908621169e-17,   -4.707367308454113e-17,   -4.517535705326342e-17,   -5.772124314976090e-17,   -8.804200044790905e-17,   -2.178465741131674e-17,   -2.881333680763293e-17;
    1.139793520796208e-16,    4.437422651548672e-18,   -3.267438741989294e-17,   -2.511165375033292e-17,   -2.838552757102909e-17,   -1.474043326635255e-17,   -1.109927918346333e-17,   -1.912739985929918e-17,   -3.715946075692235e-17,   -6.519546187633285e-17,   -4.938958652321059e-17,   -4.358765139187565e-17,   -1.049661422505236e-16,   -6.532627086844323e-17,    1.982081152310765e-17,   -1.772368330658219e-17,   -6.527211835805936e-17,   -1.346752570574594e-17,    3.058888828810693e-16;
    1.066670081255328e-16,   -7.605027718682323e-18,   -1.503078260814417e-17,    2.903493417916178e-19,   -7.266393433139062e-17,   -1.184902870265958e-17,   -3.764420416010883e-17,   -4.069249277816045e-17,   -2.099389455681444e-17,    2.145159050392920e-17,   -7.988949128970302e-17,   -4.579604924448422e-17,   -4.446687159112561e-17,   -2.039080709836938e-17,    2.643214423378448e-17,   -7.053147128843751e-17,    7.896173949070717e-17,    2.749189764728044e-17,    4.342771801890670e-19;
   -5.447178420673639e-17,    1.137631655545590e-17,   -5.500757320345621e-18,    7.584631165312438e-18,   -6.038982045008558e-17,   -4.455599300970392e-17,    1.240296792137316e-18,   -3.568858804401526e-17,   -4.265501051870825e-17,   -1.965897063194166e-17,   -2.387695415191106e-17,   -5.364362654167870e-17,   -1.106243223381823e-16,   -8.213140454196855e-17,    3.183512684696264e-17,    7.565332366642196e-18,    1.325800024778133e-17,   -8.263355277815521e-18,   -1.108445509044684e-16;
    1.997370097008705e-17,   -2.355060590986113e-17,   -5.876279551928626e-17,   -8.566286785818833e-17,   -4.873967169530952e-17,   -3.504358261907648e-18,   -7.133413326178285e-18,   -5.796150235118369e-17,   -4.496632287440822e-17,   -2.365902612710968e-17,   -1.067272491087415e-16,   -6.879571782039684e-17,   -7.972855502972470e-17,   -8.234927496852951e-17,   -2.114672640555343e-17,   -6.431969757150769e-17,   -4.200406569874332e-17,    4.246950013886419e-17,    2.462774857095104e-16;
    1.903929183093896e-16,    3.018418848199635e-19,   -6.379459474245043e-17,   -1.896391572421553e-18,    1.946299770113312e-17,    2.835730782135834e-19,   -5.588848746860068e-18,    3.511066762849900e-18,   -2.792783162391434e-17,    1.111675855536287e-17,   -4.755395093002969e-17,   -3.407583731244307e-17,   -1.183130233595262e-16,   -7.655556960931009e-17,   -7.533692975556818e-17,   -5.300135872722546e-18,   -1.423356095801246e-16,    1.482234474048383e-17,    3.518207555527116e-16;
    1.237189604236244e-16,    2.742250870824137e-17,   -2.059615498983813e-17,   -1.867169613367636e-17,   -2.231526595453115e-17,   -3.319219898934023e-17,    5.873990530091966e-18,    3.770681683556987e-17,   -4.705068799848444e-17,    2.915723218421817e-17,   -1.289051330954916e-17,   -4.995136249075216e-18,   -4.333122402555568e-18,   -9.239137233052475e-18,    3.413263595375415e-17,   -9.293173869329152e-17,   -1.178725086287136e-16,   -3.402833570476105e-17,   -1.811274654567319e-16;
   -1.245370214559508e-16,    6.349087922075108e-19,    1.532137250729280e-17,    3.054467215182072e-17,   -9.206343505466578e-18,   -4.842052323331125e-17,   -2.992138126351775e-17,   -1.088540336428162e-17,    6.419856169565495e-17,   -2.671409100873934e-17,   -6.067605321175355e-18,   -3.873010039848884e-17,   -7.252109408329745e-17,    6.912385160789956e-18,   -2.247426759125794e-17,   -5.298941895080095e-17,   -8.130239984396361e-17,    2.784231178942758e-19,    3.061060214711638e-16;
   -1.086204829423357e-16,   -2.288447209508604e-17,   -6.872964925051101e-17,   -7.795722617789935e-17,   -5.056720287725108e-17,   -4.373606511676176e-17,   -2.873271282358150e-19,   -8.841977714238447e-17,   -2.891232446598086e-17,    1.170157720720155e-17,    3.175659334022501e-17,   -2.425842729816830e-17,    1.262082750591240e-16,   -7.861615004822742e-17,   -8.362019030764417e-17,    6.835189966108989e-18,   -3.595595229975018e-17,    1.671232596756056e-17,   -6.386591830474103e-17;
    5.878093896501645e-17,    1.071018274068081e-17,   -4.903928080128052e-17,   -9.445149198351878e-18,    1.847027508695108e-17,    3.470699205462835e-17,   -1.791148406351563e-17,   -8.220866749928529e-17,   -4.716225579016499e-17,    3.005126529564972e-17,   -6.928204009299700e-17,   -1.264482361049493e-16,   -2.709435120381611e-17,    3.404340611495858e-18,   -3.495216065901343e-17,    8.321141898356282e-17,   -8.952479938471349e-17,    5.281365622611388e-18,    3.642992449316620e-16;
   -1.774442680464727e-16,   -7.938094626069870e-18,    1.888132662372643e-17,   -4.147485306582599e-17,   -2.247318000095366e-17,   -3.680129118110998e-17,   -2.728270083809936e-17,    6.325907679627374e-17,    9.448767723102546e-18,   -7.769132559509728e-17,   -8.380823162193463e-17,   -1.232549218938006e-16,    2.934387758821155e-17,   -7.343171209983623e-17,    1.250751889211865e-17,   -4.257466774959529e-17,   -1.309866386363379e-16,   -2.735311976920230e-17,   -6.807168760961102e-17;
    1.537116449206420e-17,   -2.331815296408024e-17,   -3.150128066926362e-17,    1.051747529129052e-16,    4.698170224712826e-17,   -3.122030628813222e-17,   -2.507030837811154e-17,    7.209303412494122e-17,    3.156377814824384e-17,    7.538652861682760e-17,    9.920032121592780e-17,   -1.467848466472216e-17,    7.347835651017563e-17,    2.394259920532327e-17,   -1.632857565792793e-16,    7.356613961669728e-17,   -4.931916759210765e-17,    2.983117225463516e-17,    3.043558210091748e-16;
    4.151826316555959e-16,    1.841582442096978e-17,   -9.083013475467278e-17,    9.077542320254372e-17,    1.191432884425567e-16,   -9.057208108509408e-18,    1.535938124732836e-16,   -6.236813366103378e-17,    1.235437398000236e-16,    4.519561808136174e-17,    3.232373949665218e-17,   -2.276352934274528e-17,    9.864774741432519e-17,   -4.420097100832354e-17,   -7.220474760628871e-18,   -1.578603784149757e-17,   -6.656699117705417e-17,   -5.187170137865849e-17,    3.163050062756495e-17;
   -1.387480388014787e-16,    9.572204140440023e-18,    6.793537290334349e-17,    3.634146748723172e-17,    1.368350928171360e-16,    2.583939464305407e-16,    1.756427543285390e-16,    1.139567905100378e-16,    7.333605497503691e-17,    1.648028501860521e-16,    5.529026875553643e-17,    7.453454899716133e-17,    3.821755059770990e-18,    1.762187943309124e-16,    1.476364569605228e-16,    1.401193478736337e-16,    3.931449553385390e-17,    7.392784301396560e-17,   -9.576885522132906e-17;
   -4.384981960869893e-16,   -2.664535259100378e-18,    4.789051100129171e-17,    4.290318100785839e-17,   -1.762652523940034e-17,   -1.934563620409335e-17,    3.441517903990387e-17,    3.263361803007569e-17,   -2.135965015970243e-17,   -1.654232306691483e-17,    4.375493023456300e-17,    4.850980728221543e-17,   -2.277691923957552e-18,    5.745404152435181e-18,   -3.844320695112202e-17,    8.053974154265120e-17,    3.961934946783429e-17,    4.751754545395669e-17,    1.175570057965203e-16;
   -2.332528498149612e-17,    4.391174923679131e-16,    1.920941941483453e-16,    3.424971759110815e-16,    2.148017955996493e-18,   -1.594887958679403e-16,    4.093675776779739e-16,    4.908024805799424e-16,    7.476190585391836e-16,   -4.142801553197462e-17,    7.915271865407189e-16,    3.565741316579815e-16,    8.774895381440020e-16,    2.415189901371073e-16,   -1.047169112761235e-16,    3.872320080691369e-16,    5.061124283068426e-16,    1.348824697766648e-15,    4.211162871864925e-17;
];
x{1}=[   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
   -6.250000000000000e-02,    0.000000000000000e+00,    6.250000000000000e-02,    1.250000000000000e-01,    1.875000000000000e-01,    2.500000000000000e-01,    3.125000000000000e-01,    3.750000000000000e-01,    4.375000000000000e-01,    5.000000000000000e-01,    5.625000000000000e-01,    6.250000000000000e-01,    6.875000000000000e-01,    7.500000000000000e-01,    8.125000000000000e-01,    8.750000000000000e-01,    9.375000000000000e-01,    1.000000000000000e+00,    1.062500000000000e+00;
];
x{2}=[   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02,   -6.250000000000000e-02;
    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00,    0.000000000000000e+00;
    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02,    6.250000000000000e-02;
    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01,    1.250000000000000e-01;
    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01,    1.875000000000000e-01;
    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01,    2.500000000000000e-01;
    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01,    3.125000000000000e-01;
    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01,    3.750000000000000e-01;
    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01,    4.375000000000000e-01;
    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01,    5.000000000000000e-01;
    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01,    5.625000000000000e-01;
    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01,    6.250000000000000e-01;
    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01,    6.875000000000000e-01;
    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01,    7.500000000000000e-01;
    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01,    8.125000000000000e-01;
    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01,    8.750000000000000e-01;
    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01,    9.375000000000000e-01;
    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00,    1.000000000000000e+00;
    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00,    1.062500000000000e+00;
];
