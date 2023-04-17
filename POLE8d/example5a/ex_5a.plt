 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "( 1 1 1 )" at -0.30,  1.00
set label "( 1 1 0 )" at  1.63,  1.00
set label "( 1 0 0 )" at  3.55,  1.00
set label " 1" at  4.65, -4.33
set label " 2" at  3.79, -3.43
  plot "ex_5a_CC.DAT" using 1:2 notitle  with lines lt 1 lc -1, \
 "ex_5a007.DAT" using 1:2 title "   0.7" with lines lt 1 lc  1,\
 "ex_5a010.DAT" using 1:2 title "   1.0" with lines lt 1 lc  2,\
 "ex_5a014.DAT" using 1:2 title "   1.4" with lines lt 1 lc  3,\
 "ex_5a020.DAT" using 1:2 title "   2.0" with lines lt 1 lc  4,\
 "ex_5a028.DAT" using 1:2 title "   2.8" with lines lt 1 lc  5,\
 "ex_5a040.DAT" using 1:2 title "   4.0" with lines lt 1 lc 18,\
  "ex_5a_XX.DAT" using 1:2 notitle with dots lt 0
  set term postscript eps color blacktext
  set output "ex_5a.eps"
 replot;
