 # GNUPLOT program to plot the output of POLE8.FOR
 # POLE8d developed by C. N. Tome
 set noborder
 set size ratio -1
 set key outside
 set noxtics
 set noytics
set label "( 1 1 1 )" at -0.30,  1.10
set label "( 1 1 0 )" at  1.87,  1.10
set label "( 1 0 0 )" at  4.03,  1.10
set label " 2" at  5.24, -0.00
set label " 1" at  4.27,  1.00
 plot "dotss_CC.DAT" using 1:2 notitle with lines lt 1 lc -1 ,\
 "dotss_X9.DAT" using 1:2 notitle with points   lt 1 lc -1
  set term postscript eps color blacktext
  set output "dotss.eps"
 replot;
