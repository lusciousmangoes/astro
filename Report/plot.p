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
#plot '../solver/stars.dat' using 3:($4/(3.848*10**26)) notitle lc 8 #Mass
plot '../solver/stars.dat' using 3:($4/(3.848*10**26)) notitle lc 8 #Radius
set out

