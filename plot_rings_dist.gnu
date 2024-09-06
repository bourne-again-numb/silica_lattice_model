set term postscript enhanced eps
set output 'rings_distribution.eps'
set style line 1 lt 3 lw 3 pt 3 lc rgb "red"
set style line 2 lt 3 lw 3 pt 3 lc rgb "blue"
set style line 3 lt 3 lw 3 pt 3 lc rgb "green"
set style line 4 lt 3 lw 3 pt 3 lc rgb "yellow"
set style line 5 lt 3 lw 3 pt 3 lc rgb "cyan"
set key off
set xlabel 'MC steps'
set ylabel 'rings size distribution'
plot 'rings_dist.csv' using 1:2 w l ls 2 title "3-membered rings", \
     'rings_dist.csv' using 1:3 w l ls 2 title "4-membered rings"	
