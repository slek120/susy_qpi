if (!exists("name")) name="QPI"
set title name
if (!exists("filename")) filename='test'
set output "2d/".filename.".png"

set terminal png
unset surface
unset key
set pm3d
set xlabel "{q}_{x}" 
set xrange [-pi:pi]
set ylabel "{q}_{y}"
set yrange [-pi:pi]
set view map scale 1
set size square

splot "data/".filename.".dat"

set view 130, 10, 1, 1
set pm3d scansbackward
set output "3d/".filename.".png"

splot "data/".filename.".dat"

#pause -1 "Hit return to continue"