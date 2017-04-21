set terminal gif animate delay 5

#set zr [0:50]
#set xr [0:50]
#set yr [0:50]
set xlabel 'x'
set ylabel rotate by 0 'y'
set zlabel rotate by 0 'z'
set cblabel rotate by 0 'v'
set size square
set ticslevel 0
#set cbr [0:50]
set xtics 10
set ytics 10
set ztics 10

set palette defined ( 0 'magenta', 1 'blue', 2 'cyan', 3 'green', 4 'yellow', 5 'orange', 6 'red' ) 

#set output 'Universe3D.gif'
#do for [i=0:9] { splot sprintf('../ControlData/universe0000%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}
#do for [i=10:99] { splot sprintf('../ControlData/universe000%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}
#do for [i=100:392] { splot sprintf('../ControlData/universe00%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}
#set out

set terminal epslatex size 5.33in,4.0in

set xlabel '$x$'
set ylabel rotate by 0 '$y$'
set zlabel rotate by 0 '$z$'
set cblabel rotate by 0 '$v$'

#set output 'Control.tex'
#plot '../Results/Control.dat' using 2:3:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
#set out

set cbr [0:50]
#set output 'N64.tex'
#splot '../Results/N64.dat' using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
#set out

set cbr [0:100]
set output '3mass.tex'
splot '../Results/m3l0k0.dat' using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
set out

set cbr [0:10]
#set output '1Ncells.tex'
#splot '../Results/1Ncells.dat' using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
#set out

set cbr [0:7.5]
set output 'Ordered.tex'
splot '../Results/Ordered.dat' using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
set out

set xtics 0.2
set ytics 0.2
set ztics 0.2
set cbr [0:7.5]
set output '1boxsize.tex'
splot '../Results/1boxsize.dat' using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle
set out


#set ticslevel 0
unset tics
unset border
unset xlabel
unset ylabel
unset cblabel
unset colorbox

#set output 'Rotate.gif'

speed=2

#do for [i=0:9] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('../ControlData/universe0000%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}
#do for [i=10:99] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('../ControlData/universe000%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}
#do for [i=100:392] { set view 60, (speed*i-floor(speed*i/360)*360); splot sprintf('../ControlData/universe00%d.dat', i) using 2:3:4:(sqrt($5*2 + $6**2 + $7**2)) with dots palette notitle}

#set out
