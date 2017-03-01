set terminal pngcairo size 700,524 enhanced font 'Verdana,10'

set notitle
set xlabel 'T (K)'
set ylabel 'L (Solar Lumens)'

set logscale xyz
set format x "10^%T"
set format y "10^%T"

set output sprintf('./HR.png')
set grid front
plot '../solver/stars.dat' using 3:($4/(3.848*10**26)):($2/10**26) with circles notitle #Mass
#plot '../solver/stars.dat' using 3:($4/(3.848*10**26)):($1/10**4) with circles notitle #Radius
set out
unset grid

