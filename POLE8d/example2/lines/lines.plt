 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "(001)" at  0.00,  1.04
set label "(001)" at  3.05, -0.10
set label "(110)" at  5.13, -0.10
set label "(111)" at  5.13,  1.98
  plot "lines_CC.DAT" using 1:2 notitle  with lines lt 1 lc -1, \
 "lines005.DAT" using 1:2 title "   0.5" with lines lt 1 lc  1,\
 "lines010.DAT" using 1:2 title "   1.0" with lines lt 1 lc  2,\
 "lines020.DAT" using 1:2 title "   2.0" with lines lt 1 lc  3,\
 "lines040.DAT" using 1:2 title "   4.0" with lines lt 1 lc  4,\
 "lines080.DAT" using 1:2 title "   8.0" with lines lt 1 lc  5,\
 "lines160.DAT" using 1:2 title "  16.0" with lines lt 1 lc 18,\
 "lines320.DAT" using 1:2 title "  32.0" with lines lt 1 lc  7,\
  "lines_XX.DAT" using 1:2 notitle with dots lt 0
  set term postscript eps color blacktext
  set output "lines.eps"
 replot;
