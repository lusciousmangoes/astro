#set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set terminal epslatex color size 6.0in,4.5in

set notitle
set xlabel '$T$ (K)'
set ylabel '$L$ ($L_\odot$)'

set logscale xyz
set format x '$10^{%T}$'
set format y '$10^{%T}$'


#set logscale cb
#set cbrange[10**6:10**9]

set xr [*:*] reverse
set key below #bottom left

set output sprintf('./badHR.tex')
set grid front
plot '../solver/stars_linear.dat' using 3:($4/(3.848*10**26)):2 lc 8 lt 7 lw 1 notitle, \
'../solver/stars_log.dat' using 3:($4/(3.848*10**26)):2 lc 8 lt 7 lw 1 notitle
set out

