#set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set terminal epslatex color size 6.0in,4.5in

set notitle
set xlabel '$T$ (K)'
set ylabel '$L$ ($L_\odot$)'

set logscale xyz
set format x '$10^{%T}$'
set format y '$10^{%T}$'

set xr [*:*] reverse
set key below #bottom left

set output sprintf('./HR.tex')
set grid front
plot '../MESA-Web_job_1-solar-mass-2percent-metalicity/HR.dat' using 2:1 notitle with linespoints lc 'cyan' lt -1 lw 5, \
'../MESA-Web_Job_0226165219/HR.dat' using 2:1 notitle with linespoints lc 'green' lt -1 lw 5, \
'../MESA-Web_job_point4-solar-mass-point4percent-metalicity/HR.dat' using 2:1 notitle with linespoints lc 'yellow' lt -1 lw 5, \
'../MESA/HR.dat' using 2:1 notitle with linespoints lc 'orange' lt -1 lw 5, \
'../MESA2/HR.dat' using 2:1 notitle with linespoints lc 'red' lt -1 lw 5, \
'../MESA3/HR.dat' using 2:1 notitle with linespoints lc 'navy' lt -1 lw 5, \
'../MESA4/HR.dat' using 2:1 title 'MESA star' with linespoints lt -1 lw 5, \
'../MESA5/HR.dat' using 2:1 notitle with linespoints lc 'brown' lt -1 lw 5, \
'../solver/real_stars.dat' using 1:2 title 'Real Stars' lc 'magenta' lt 7 lw 1 ps 1.7, \
'../solver/stars.dat' using 3:($4/(3.848*10**26)) title 'Simulated Stars' lc 8 lt 7 lw 1
set out

