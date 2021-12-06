set size 0.95, 0.95
set xrange [95:135]  
set yrange [40:60]  
set bmargin 5
set lmargin 12

set title "LaBr/CeBr resolution (at RoSPHERE)"

set xtics font ",14"
set ytics font ",14"  

set encoding iso_8859_1 
set xlabel "Domain" font ",16" offset 0,-1
# set xlabel " "
# set ylabel " "
set ylabel "Resolution, keV" font ",16" offset -2,0

set style line 1 linecolor rgb '#ff0000' linetype 2 linewidth 3 
set style line 110 linecolor '#black' linetype 2 linewidth 2 dt 3 pointtype 7 pointsize 1.5 

set style line 112 linecolor '#black' linetype 1 linewidth 3 
set style line 114 linecolor rgb '#ff8040' linetype 2 linewidth 1 dt 3 pointtype 9 pointsize 2.5 
#set logscale x
#set xtics format " " 
# set xtics (10, 100, 1000)
set xtics autofreq

#set label 1 at 13, 15
#set label 1 '(b)' font ",18" tc '#black'

set key font ",14" at 135.0, 59.
#set nokey
plot 'resolution_1173.txt'  using 1:6 linestyle 110  linecolor rgb '#0000FF' title "1173 keV",\
     'resolution_1332.txt'  using 1:6 linestyle 110  linecolor rgb '#8A2BE2' title "1332 keV" \
#with lines, \
#     'intensity.dat' using 1:2:($3-$2)  linestyle 110 linecolor '#black' linetype 2 notitle  with yerrorbars, \
#     'intensity.dat' using 1:2    linestyle 110 notitle,  
#     'tau15.dat' using 4:6    linestyle 112 linewidth 1 notitle with lines, \
#     'tau15.dat' using 4:7    linestyle 112 linewidth 1 notitle with lines \
#     'e2plot.dat' using 1:3    linestyle 112  notitle with linespoints, \
#     'e2plot.dat' using 1:3    linestyle 112  title "^{112}Te", \
#     'e2plot.dat' using 1:5:6  linestyle 114  linetype 2 linecolor rgb '#ff8040' notitle  with yerrorbars, \
#     'e2plot.dat' using 1:5    linestyle 114  notitle  with linespoints, \
#     'e2plot.dat' using 1:5    linestyle 114  title "^{114}Te", \
#     'e2plot.dat' using 1:2    linestyle 1 title "TAC"  with line

