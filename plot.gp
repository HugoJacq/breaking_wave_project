# Plot the heatmap
set pm3d map
set view map
set xlabel "Time (s)"
set ylabel "Layer"
set cblabel "Value"
set yrange [0:15]
set terminal pngcairo size 800,600 enhanced font 'Verdana,12'
set output 'u_profile.png'
set size 0.9, 0.9
splot "u_profile.dat" using 1:2:3 with pm3d
unset output

