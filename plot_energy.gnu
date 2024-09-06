set term postscript enhanced eps
set output 'energy.eps'
set style line 1 lt 4 lw 1 pt 3 lc rgb "black"
set key off
set xlabel 'MC steps'
set ylabel 'Energy per lattice site'
plot 'lattice_energy.csv' using 1:3 w lp ls 1

