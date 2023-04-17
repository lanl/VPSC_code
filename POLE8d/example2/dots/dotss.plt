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
 plot "dotss_CC.DAT" using 1:2 notitle with lines lt 1 lc -1 ,\
 "dotss_X9.DAT" using 1:2 notitle with points   lt 1 lc -1
  set term postscript eps color blacktext
  set output "dotss.eps"
 replot;
