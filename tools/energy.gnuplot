#!/usr/local/bin/gnuplot
set terminal pdf color solid font "helvetica,20"
set output "energy.pdf"
set xlabel "Simulation Steps"
set ylabel "Simulated Energy(eV/atom)"
set xtics rotate by 45 right
plot 'kmc_log.txt' u 1:($3/10976) w l lw 4 notitle
