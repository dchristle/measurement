set style data lines
set x2label ""
set xlabel "X"
set y2label ""
set ylabel "Y"
set grid
plot "202606_127327641.tmp" using ($1+0.000000+0.000000*column(-1)) every ::0 binary format='%float64' title "tmp/202606_127327641.tmp, Y vs " axes x1y1