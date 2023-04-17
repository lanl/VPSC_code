 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "( 0 0 0 1 )" at -0.30,  1.10
set label "( 1 0-1 1 )" at  1.87,  1.10
set label "( 1 0-1 0 )" at  4.03,  1.10
set label " 1" at  5.24, -2.44
set label " 2" at  4.27, -1.43
  plot "exam3_CC.DAT" using 1:2 notitle  with lines lt 1 lc -1, \
 "exam3007.DAT" using 1:2 title "   0.7" with lines lt 1 lc  1,\
 "exam3010.DAT" using 1:2 title "   1.0" with lines lt 1 lc  2,\
 "exam3014.DAT" using 1:2 title "   1.4" with lines lt 1 lc  3,\
 "exam3020.DAT" using 1:2 title "   2.0" with lines lt 1 lc  4,\
 "exam3028.DAT" using 1:2 title "   2.8" with lines lt 1 lc  5,\
 "exam3040.DAT" using 1:2 title "   4.0" with lines lt 1 lc 18,\
  "exam3_XX.DAT" using 1:2 notitle with dots lt 0
  set term postscript eps color blacktext
  set output "exam3.eps"
 replot;
