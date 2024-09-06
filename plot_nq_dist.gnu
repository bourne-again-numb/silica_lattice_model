set term postscript enhanced eps
set key on

set style line 1 lt 5 lw 2 pt 3 lc rgb "red"
set style line 2 lt 5 lw 2 pt 3 lc rgb "blue"
set style line 3 lt 5 lw 2 pt 3 lc rgb "green"
set style line 4 lt 5 lw 2 pt 3 lc rgb "magenta"
set style line 5 lt 5 lw 2 pt 3 lc rgb "cyan"
set style line 6 lt 5 lw 2 pt 3 lc rgb "black"

set output 'nq_mcsteps.eps'
set logscale x
set xlabel 'MC steps'
set ylabel 'Relative Qn distribution'
set xrange[20:1e+7]
set yrange[0:1]
plot 'qn.c_mcsteps.csv' using 1:2 w l ls 1 title 'Q0',\
     'qn.c_mcsteps.csv' using 1:3 w l ls 2 title 'Q1',\
     'qn.c_mcsteps.csv' using 1:4 w l ls 3 title 'Q2',\
     'qn.c_mcsteps.csv' using 1:5 w l ls 4 title 'Q3',\
     'qn.c_mcsteps.csv' using 1:6 w l ls 5 title 'Q4',\
     'qn.c_mcsteps.csv' using 1:7 w l ls 6 title 'c'
	
set output 'nq_c.eps'
unset logscale x
set xlabel 'Degree of Condensation'
set ylabel 'Relative Qn distribution'
set xrange[0:1]
set yrange[0:1]
plot 'qn.c_mcsteps.csv' using 7:2 w l ls 1 title 'Q0',\
     'qn.c_mcsteps.csv' using 7:3 w l ls 2 title 'Q1',\
     'qn.c_mcsteps.csv' using 7:4 w l ls 3 title 'Q2',\
     'qn.c_mcsteps.csv' using 7:5 w l ls 4 title 'Q3',\
     'qn.c_mcsteps.csv' using 7:6 w l ls 5 title 'Q4'

