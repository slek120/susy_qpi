set terminal png
set output 'test.png'
#set title "Test" 
set xlabel "qx" 
set ylabel "qy"
set view map scale 1
set xrange [-pi:pi]
set yrange [-pi:pi]
set size square
unset surface 
set pm3d

splot "< sqlite3 -column data.db 'select distinct qx,qy,absresult from susy_qpi order by qx, qy;' | awk -f add_blanks.awk"
#pause -1 "Hit return to continue"