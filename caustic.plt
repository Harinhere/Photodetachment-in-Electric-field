reset
set key top right font ",12"
set multiplot layout 2,1
set lmargin at screen 0.10; set rmargin at screen 0.90
set tmargin at screen 0.90; set bmargin at screen 0.55

set ylabel "z (a.u.)" font ",12"
unset xtics
set label "Classically\n forbidden\n region" at graph 0.02,0.4 font ",12"
set yr[10:-3]
set xr[-30:30]
set ytics font ",12"
plot for [i=101:181] 'stat_cstcs'.i.'.out' u 2:3 w l lc -1  lw 0.5 title "",\
'caustic.out' u 1:2 w l lc 7 lw 2 title "Caustic" 

set lmargin at screen 0.10; set rmargin at screen 0.90
set tmargin at screen 0.55; set bmargin at screen 0.20
unset yr
set xtics
set xtics font ",12"
unset ylabel
set ylabel "Re [G]" font ",12"
set xlabel "x (a.u.)" font",12"
unset label

set yr[-0.05:0.05]

p 'green_el.out' u 1:2 w l lc 7 lw 2 dt 2 title "QM",\
'green_el.out' u 1:4 w l lc -1 lw 1.5  title "SC"

unset multiplot