#set terminal qt
set terminal png
set output 'test.png'
#set title "Test" 
#set xlabel "qx" 
#set ylabel "qy" 
set view map scale 1
set xrange [0:pi]
set yrange [0:pi]
unset surface 
set pm3d

splot 'w=.00+1.00i_susy_qpi.dat'
pause -1 "Hit return to continue"