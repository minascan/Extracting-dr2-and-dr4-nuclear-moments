#!/bin/sh
#Trans.1 : 1s(2).2s_2S.2p_1P 1-  --> 1s(2).2p(2)3P2_3P 0+   lin_par = -6.0413444
#Trans.2 : 1s(2).2p(2)1S0_1S 0+  --> 1s(2).2s_2S.2p_3P 1-   lin_par =  7.7281533

./extract_dr2_dr4 <<EOF
92   # atomic number
238  # mass number of isotope A
236  # mass number of isotope A'
true # logic variable controling the heavier isotope
1s(2).2p(2)_3P  #Trans.1
5708068.97
5709843.65
1s(2).2s2p_1P
36358007.94
36359001.82
1s(2).2s.2p_3P  #Trans.2
2401822.81
2402710.12
1s(2).2p(2)_1S
72727658.08
72729668.71
5723257.74  #Trans.1 guessed rms radius 236
36366513.08
2409416.24  #Trans.2 guessed rms radius 236
72744864.63
5.8571  #Values of rms radii
5.8431
5.7363
1.5495363D+05  #Electronic factors from RIS4 Trans.1
-2.0435053D+02
5.3676241D-01
-8.8239448D-04
-2.2251995D+05  #Electronic factors from RIS4 Trans.2
2.8944623D+02
-7.6046445D-01
1.2468882D-03
0.001  #factor to be multiplied with the pseudo-FS data
EOF
