#set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set terminal epslatex color size 6.0in,4.5in

set notitle
set xlabel '$T$ (K)'
set ylabel '$L$ ($L_\odot$)'

set logscale xyz
set format x '$10^{%T}$'
set format y '$10^{%T}$'

set xr [*:*] reverse

set output sprintf('./HR.tex')
set grid front
plot '../MESA-Web_Job/HR.dat' using 2:1 title 'MESA star' lc 'green' lt 7 lw 1, \
'../MESA-Web_job_1-solar-mass-2percent-metalicity/HR.dat' using 2:1 title 'MESA star' lc 'blue' lt 7 lw 1, \
'../MESA-Web_Job_100_solar_mass/HR.dat' using 2:1 title 'MESA star' lc 'yellow' lt 7 lw 1, \
'../MESA-Web_Job_0226165219/HR.dat' using 2:1 title 'MESA star' lc 'cyan' lt 7 lw 1, \
'../MESA-Web_job_point4-solar-mass-point4percent-metalicity/HR.dat' using 2:1 title 'MESA star' lc 'purple' lt 7 lw 1, \
'../solver/real_stars.dat' using 1:2 title 'Real Stars' lc 'red' lt 7 lw 1, \
'../solver/stars.dat' using 3:($4/(3.848*10**26)) title 'Simulated Stars' lc 8 lt 7 lw 1
set out

