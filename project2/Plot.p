set terminal gif animate delay 5

set zr [0:50]
set xr [0:50]
set yr [0:50]
set xlabel 'x'
set ylabel 'y'
set size square

#set output 'Universe3D_take2.gif'
#do for [i=0:9] { splot sprintf('Data/universe0000%d.dat', i) using 2:3:4 with dots notitle}
#do for [i=10:99] { splot sprintf('Data/universe000%d.dat', i) using 2:3:4 with dots notitle}
#do for [i=100:392] { splot sprintf('Data/universe00%d.dat', i) using 2:3:4 with dots notitle}
#set out

set output 'UniverseJoe_take2.gif'
do for [i=0:9] { plot sprintf('Data/universe0000%d.dat', i) using 2:3 with dots notitle}
do for [i=10:99] { plot sprintf('Data/universe000%d.dat', i) using 2:3 with dots notitle}
do for [i=100:980] { plot sprintf('Data/universe00%d.dat', i) using 2:3 with dots notitle}
set out
