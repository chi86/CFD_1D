set encoding utf8

# See https://github.com/Gnuplotting/gnuplot-palettes
# Line styles (colorbrewer Set1)
set style line 1 lc rgb '#E41A1C' pt 1 ps 1 lt 1 lw 4 # red
set style line 2 lc rgb '#377EB8' pt 2 ps 1 lt 1 lw 4 # blue
set style line 3 lc rgb '#4DAF4A' pt 6 ps 1 lt 1 lw 4 # green
set style line 4 lc rgb '#984EA3' pt 3 ps 1 lt 1 lw 4 # purple
set style line 5 lc rgb '#FF7F00' pt 4 ps 1 lt 1 lw 4 # orange
set style line 6 lc rgb '#FFFF33' pt 5 ps 1 lt 1 lw 4 # yellow
set style line 7 lc rgb '#A65628' pt 7 ps 1 lt 1 lw 4 # brown
set style line 8 lc rgb '#F781BF' pt 8 ps 1 lt 1 lw 4 # pink
# Palette
set palette maxcolors 8
set palette defined ( 0 '#E41A1C', 1 '#377EB8', 2 '#4DAF4A', 3 '#984EA3',\
4 '#FF7F00', 5 '#FFFF33', 6 '#A65628', 7 '#F781BF' )

# Standard border
set style line 11 lc rgb '#808080' lt 1 lw 3
set border 0 back ls 11
set tics out nomirror

# Standard grid
set style line 12 lc rgb '#808080' lt 0 lw 1
set grid back ls 12
unset grid


ReTau=100

### laminar channel
f(x) = ReTau/1*(0.25-x*x)
### laminar pipe
# f(x) = ReTau/1*(0.25-x*x)



set terminal pdfcairo enhanced color dashed font "Alegreya, 14" \
rounded size 16 cm, 9.6 cm
set output 'laminarFlow.pdf'

# Default encoding, line styles, pallette, border and grid are set in
# /usr/local/share/gnuplot/x.y/gnuplotrc.

set xlabel "y^+"
set ylabel "L/T"
set grid
set key left top
#set xrange[0:0.5]
#set yrange[-1:1]

set logscale x
#set logscale y

plot "data.dat" u 4:5 w l ls 1 title "data", \
     "data.dat" u 4:(f($2)) every 10 w p ls 2 title "laminar flow"