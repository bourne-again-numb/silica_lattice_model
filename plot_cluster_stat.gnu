set term postscript enhanced eps
set output 'cluster_max.eps'
set style line 1 lt 1 lw 3 pt 3 lc rgb "black"
set key off
set xrange[10:10e+06]
#set yrange[0:100]
set xlabel 'MC steps'
set ylabel 'Maximumm Cluster Size'
plot 'cluster_stat.csv' using 1:2 w l ls 1

set output 'cluster_avg.eps'
set style line 1 lt 1 lw 3 pt 3 lc rgb "black"
set key off
set xrange[10:10e+06]
#set yrange[0:100]
set xlabel 'MC steps'
set ylabel 'Average Cluster Size'
plot 'cluster_stat.csv' using 1:3 w lp ls 1

set output 'cluster_num.eps'
set style line 1 lt 1 lw 3 pt 3 lc rgb "black"
set key off
set xrange[10:10e+06]
#set yrange[0:100]
set xlabel 'MC steps'
set ylabel 'Number of Clusters'
plot 'cluster_stat.csv' using 1:4 w lp ls 1

