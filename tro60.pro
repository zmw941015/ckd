c     tropical  :   p(mb)  t(k)  h2o(g/g)  o3(g/g)
c     np=60, surface temp = 300.0
c
      data (plevel(i),i=1,np+1)/
     &    0.0712906,   0.1000000,   0.1402710,   0.1967600,   0.2759970,
     &    0.3871,    0.5431,    0.7617,    1.0685,    1.4988,    2.1024,
     &    2.9490,    4.1366,    5.8025,    8.1392,   11.4170,   16.0147,
     &   22.4640,   31.5105,   44.2001,   62.0000,   85.7750,  109.5500,
     &  133.3250,  157.1000,  180.8750,  204.6500,  228.4250,  252.2000,
     &  275.9750,  299.7500,  323.5250,  347.3000,  371.0750,  394.8500,
     &  418.6250,  442.4000,  466.1750,  489.9500,  513.7250,  537.5000,
     &  561.2750,  585.0500,  608.8250,  632.6000,  656.3750,  680.1500,
     &  703.9250,  727.7000,  751.4750,  775.2500,  799.0250,  822.8000,
     &  846.5750,  870.3500,  894.1250,  917.9000,  941.6750,  965.4500,
     &  989.2250, 1013.0000/
c
      data (tlayer(i),i=1,np)/
     &  225.8594,  232.3394,  239.0249,
     &  245.8975,  252.9619,  260.2183,  266.9444,  269.2215,  267.0145,
     &  263.0619,  257.6949,  252.0602,  246.5490,  241.1625,  235.8975,
     &  230.7513,  225.7210,  220.8037,  215.8554,  209.3228,  201.4490,
     &  196.9717,  200.1833,  207.3744,  213.7087,  219.3563,  224.4651,
     &  229.1386,  233.4520,  237.4624,  241.2137,  244.7407,  248.0714,
     &  251.2287,  254.2318,  257.0965,  259.8364,  262.4629,  264.9861,
     &  267.4148,  269.7563,  272.0171,  274.2024,  276.3143,  278.3489,
     &  280.2858,  282.0777,  283.6658,  285.0422,  286.2809,  287.4879,
     &  288.7410,  290.0624,  291.4273,  292.8017,  294.1639,  295.5044,
     &  296.8197,  298.1095,  299.3740/
c
      data (wlayer(i),i=1,np)/
     &  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,
     &  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,
     &  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,
     &  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,  0.3250E-05,
     &  0.3250E-05,  0.3250E-05,  0.3254E-05,  0.3379E-05,  0.3885E-05,
     &  0.5967E-05,  0.1417E-04,  0.2944E-04,  0.5517E-04,  0.9524E-04,
     &  0.1538E-03,  0.2351E-03,  0.3388E-03,  0.4425E-03,  0.5666E-03,
     &  0.7130E-03,  0.8837E-03,  0.1081E-02,  0.1305E-02,  0.1560E-02,
     &  0.1846E-02,  0.2165E-02,  0.2520E-02,  0.2910E-02,  0.3339E-02,
     &  0.3807E-02,  0.4314E-02,  0.4882E-02,  0.5916E-02,  0.7010E-02,
     &  0.8139E-02,  0.9275E-02,  0.9960E-02,  0.1054E-01,  0.1116E-01,
     &  0.1184E-01,  0.1258E-01,  0.1338E-01,  0.1423E-01,  0.1516E-01/
c
      data (olayer(i),i=1,np)/
     &  0.1288E-05,  0.1622E-05,  0.1975E-05,  0.2347E-05,  0.2739E-05,
     &  0.3151E-05,  0.3574E-05,  0.4269E-05,  0.5379E-05,  0.6905E-05,
     &  0.8719E-05,  0.9931E-05,  0.1045E-04,  0.1122E-04,  0.1234E-04,
     &  0.1212E-04,  0.1028E-04,  0.7945E-05,  0.5155E-05,  0.2594E-05,
     &  0.1040E-05,  0.3939E-06,  0.2292E-06,  0.1898E-06,  0.1641E-06,
     &  0.1454E-06,  0.1282E-06,  0.1139E-06,  0.1023E-06,  0.9333E-07,
     &  0.8689E-07,  0.8190E-07,  0.7752E-07,  0.7416E-07,  0.7200E-07,
     &  0.7025E-07,  0.6863E-07,  0.6710E-07,  0.6568E-07,  0.6435E-07,
     &  0.6307E-07,  0.6188E-07,  0.6077E-07,  0.5968E-07,  0.5899E-07,
     &  0.5869E-07,  0.5839E-07,  0.5798E-07,  0.5737E-07,  0.5669E-07,
     &  0.5602E-07,  0.5531E-07,  0.5452E-07,  0.5374E-07,  0.5297E-07,
     &  0.5213E-07,  0.5115E-07,  0.5012E-07,  0.4910E-07,  0.4810E-07/

