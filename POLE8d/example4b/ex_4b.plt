 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "( 0 0 1 )" at -0.30,  0.65
set label "( 1 0 0 )" at  0.78,  0.65
set label "( 1 1 0 )" at  1.87,  0.65
set label "( 1 0 1 )" at  2.95,  0.65
set label "( 1 0 2 )" at  4.03,  0.65
set label "( 1 1 2 )" at  5.12,  0.65
set label " 1" at  5.87, -2.44
set label " 2" at  5.36, -1.89
  plot "ex_4b_CC.DAT" using 1:2 notitle  with lines lt 1 lc -1, \
 "ex_4b007.DAT" using 1:2 title "   0.7" with lines lt 1 lc  1,\
 "ex_4b010.DAT" using 1:2 title "   1.0" with lines lt 1 lc  2,\
 "ex_4b014.DAT" using 1:2 title "   1.4" with lines lt 1 lc  3,\
 "ex_4b020.DAT" using 1:2 title "   2.0" with lines lt 1 lc  4,\
 "ex_4b028.DAT" using 1:2 title "   2.8" with lines lt 1 lc  5,\
 "ex_4b040.DAT" using 1:2 title "   4.0" with lines lt 1 lc 18,\
 "ex_4b056.DAT" using 1:2 title "   5.7" with lines lt 1 lc  7,\
  "ex_4b_XX.DAT" using 1:2 notitle with dots lt 0
  set term postscript eps color blacktext
  set output "ex_4b.eps"
 replot;
