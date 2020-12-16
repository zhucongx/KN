#!/usr/local/bin/gnuplot
set terminal pdf color solid font "helvetica,20"
set output "time.pdf"
set xlabel "Simulation Steps"
set ylabel "Simulated Time(s)"
set xtics rotate by 45 right

plot 'kmc_log.txt' u 1:2 w l lw 4 notitle
