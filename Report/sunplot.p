#set terminal pngcairo size 700,524 enhanced font 'Verdana,10'
set terminal epslatex color size 6.2in,2.7in

set notitle
set xlabel '$r$ ($R_{\odot}$)'
set grid front

set output sprintf('./sunM.tex')
set ylabel '$M$ ($M_\odot$)'
plot '../solver/sun.dat' using ($1/(6.957*10**8)):($2/(2.0*10**30)) notitle with lines lc 8
set out

set output sprintf('./sunT.tex')
set ylabel '$T$ (K)'
plot '../solver/sun.dat' using ($1/(6.957*10**8)):3 notitle with lines lc 8
set out

set output sprintf('./sunL.tex')
set ylabel '$L$ ($L_\odot$)'
plot '../solver/sun.dat' using ($1/(6.957*10**8)):($4/(3.848*10**26)) notitle with lines lc 8
set out

