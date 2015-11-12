plot 'DatenniceT2.txt' matrix with image
set palette grey
set xrange [0:200]
set yrange [0:200]
set term png            
set output "spinlatticeT2.png" 
replot
set term x11
