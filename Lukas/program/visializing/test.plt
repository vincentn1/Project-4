plot 'Daten.txt' matrix with image
set palette grey
set xrange [0:100]
set yrange [0:100]
set term png            
set output "spinlattice.png" 
replot
set term x11
