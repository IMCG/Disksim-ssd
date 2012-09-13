set terminal postscript 
set output "page.ps"
set noclip points
set clip one
set noclip two
set border
set boxwidth
set dummy x,y
set format x "%g"
set format y "%g"
set format z "%g"
set grid
set key below #9,60
set nolabel
set noarrow
set nologscale
set offsets 0, 0, 0, 0
set nopolar
set angles radians
set noparametric
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
set nocontour
set clabel
set nohidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size 1,1
set data style points
set function style lines
set xzeroaxis
set yzeroaxis
set tics in
set ticslevel 0.5
set xtics
set ytics
set ztics
set title "Execution times" 0,0
set notime
set rrange [-0 : 10]
set trange [-5 : 5]
set urange [-5 : 5]
set vrange [-5 : 5]
set xlabel "Barnes" 0,0
set xrange [1:7]
set ylabel "cycles" 0,0
set yrange [1 : 2.5]
set zlabel "" 0,0
set zrange [-10 : 10]
set autoscale r
set autoscale t
set autoscale xy
set autoscale z
set zero 1e-08
set pointsize 2
#plot "geneid.quad_pro.data" using 1:2 with imp linewidth 40
#set xtics ("1" 1, "2" 2, "3" 3, "4" 4, "5" 5, "6" 6, "7" 7, "8" 8, "9" 9, "10" 10, "11" 11, "12" 12, "13" 13, "14" 14, "15" 15, "16" 16, "17", 17, "18" 18, "19"19, "20" 20)
#set xtics ("1" 1, "2" 2, "3" 3, "4" 4, "5" 5, "6" 6, "7" 7, "8" 8, "9" 9)
#plot "barnes" with boxes fill solid 0.00
plot "page" with boxes fill solid 0.00
#plot "1" with boxes fill solid 0.00, "2" with boxes fill solid 0.10, "3" with boxes fill solid 0.20, "4" with boxes fill solid 0.30, "5" with boxes fill solid 0.40, "6" with boxes fill solid 0.50
