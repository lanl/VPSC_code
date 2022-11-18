 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "( 0 0 0 1 )" at -0.30,  0.88
set label "( 1 0-1 0 )" at  1.32,  0.88
set label "( 2-1-1 0 )" at  2.95,  0.88
set label "( 1 0-1 1 )" at  4.57,  0.88
set label " 1" at  5.55, -3.66
set label " 3" at  4.82, -2.88
  plot "ex_5b_CC.DAT" using 1:2 notitle  with lines lt 1 lc -1, \
 "ex_5b007.DAT" using 1:2 title "   0.7" with lines lt 1 lc  1,\
 "ex_5b010.DAT" using 1:2 title "   1.0" with lines lt 1 lc  2,\
 "ex_5b014.DAT" using 1:2 title "   1.4" with lines lt 1 lc  3,\
 "ex_5b020.DAT" using 1:2 title "   2.0" with lines lt 1 lc  4,\
 "ex_5b028.DAT" using 1:2 title "   2.8" with lines lt 1 lc  5,\
  "ex_5b_XX.DAT" using 1:2 notitle with dots lt 0
  set term postscript eps color blacktext
  set output "ex_5b.eps"
 replot;
