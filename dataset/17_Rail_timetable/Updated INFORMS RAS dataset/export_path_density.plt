set title "Dynamic Density Contour (Path 1 ) Unit: veh/mile/lane" 
set xlabel "Time Horizon"
set ylabel "Space (Node Sequence)"  offset -1
set xtics (" 0:00" 0 ," 1:00" 60 ," 2:00" 120 ," 3:00" 180 ," 4:00" 240 ," 5:00" 300 ," 6:00" 360 ," 7:00" 420 ," 8:00" 480 ," 9:00" 540 ) 
set ytics ("null" 0, "null" 180, "null" 190, "null" 200, "null" 400, "null" 580, "null" 590, " " 600)
set xrange [0:541] 
set yrange [0:600] 
set palette defined (0 "white", 10 "green", 30 "yellow", 50 "red")
set pm3d map
splot 'C:\GitHub\learning-transportation-engineering-and-traffic-analysis\undergraduate_student_project\data_sets_GMNS0.9\17_Rail_INFORMS_RAS\Updated ras dataset\export_path_density.txt' matrix notitle
