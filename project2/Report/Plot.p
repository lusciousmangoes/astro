set terminal gif animate delay 5

set zr [0:50]
set xr [0:50]
set yr [0:50]
set xlabel '$x$'
set ylabel rotate by 0 '$y$'
set cblabel rotate by 0 '$z$'
set size square

set palette defined ( 0 'magenta', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 'orange', 6 'red' ) 

#set output 'Universe2D.gif'
#do for [i=0:9] { plot sprintf('Data/universe0000%d.dat', i) using 2:3:4 with dots palette notitle}
#do for [i=10:99] { plot sprintf('Data/universe000%d.dat', i) using 2:3:4 with dots palette notitle}
#do for [i=100:392] { plot sprintf('Data/universe00%d.dat', i) using 2:3:4 with dots palette notitle}
#set out

set terminal epslatex size 5.33in,4.0in

set ticslevel 0

set output 'Control.tex'
plot '../Results/Control.dat' using 2:3:4 with dots palette notitle
set out

set output 'Control3D.tex'
splot '../Results/Control.dat' using 2:3:4 with dots palette notitle
set out


#set ticslevel 0
#unset tics
#unset border
#unset xlabel
#unset ylabel
#unset cblabel
#unset colorbox

#set output 'Rotate.gif'

#speed=2

#do for [i=0:9] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('Data/universe0000%d.dat', i) using 2:3:4 with dots palette notitle}
#do for [i=10:99] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('Data/universe000%d.dat', i) using 2:3:4 with dots palette notitle}
#do for [i=100:392] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('Data/universe00%d.dat', i) using 2:3:4 with dots palette notitle}

#set out
