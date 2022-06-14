# Simple script to make boxplots
# Assumes column headers are group labels
# run gnuplot -e "filename='/path/to/file.txt'" boxplot.gnuplot


set xtics auto
set yrange [*:*]
set style boxplot outliers pointtype 7
set style fill solid 0.25 border -1
set style data boxplot
set boxwidth 0.1
set pointsize 0.5


factors = system('head -1 '.filename)
n_f = words(factors)
set xtic ("" 1)
set for [i=1:n_f] xtics add (word(factors,i) i)


plot for [i=1:n_f] filename using (i):i
pause -1
